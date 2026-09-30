%% Image analysis of bacterial accumulation around spheroids from time-lapse sequences
%
% This script quantifies the spatial and temporal dynamics of bacterial
% accumulation around spheroids from time-lapse microscopy images. The
% spheroids were co-inoculated with two differentially labelled bacterial
% strains.
%
% Workflow:
% 1) Identify the spheroid boundary in the bright-field channel.
% 2) Quantify fluorescence intensity around and inside the spheroid
%    in two fluorescent channels.
% 3) Generate radial fluorescence intensity profiles over time
%    for the two fluorescent channels.
% 4) Identify angular sectors with pronounced dispersal waves
%    in the two fluorescent channels.
%
% The spheroid position is identified independently in each frame,
% and the spheroid center is estimated by fitting a circle to its boundary.

clear all
close all
clc

% microscopy details
pixelsize = 2.77; % in microns
timeframeint = 30/60; % time between frames in minutes

% load the bright-field movie
filenameBF = 'example_BF.tif'; % CHANGE FILE NAME (bright-field channel)
% load the fluorescent movies
filenameGFP = 'example_GFP.tif'; % CHANGE FILE NAME (fluorescent channel N1)
filenamemCh = 'example_mCh.tif'; % CHANGE FILE NAME (fluorescent channel N2)

infoBF = imfinfo(filenameBF);      % read metadata
numFrames = numel(infoBF);         % number of frames

movieDataBF = zeros(infoBF(1).Height, infoBF(1).Width, numFrames, 'uint8');

for k = 1:numFrames
    movieDataBF(:,:,k) = imread(filenameBF, k);
end

%% Identify the position of the spheroid for each frame

boundaries = cell(numFrames,1);
beltMasks = cell(numFrames,1); % surface
insideMasks = cell(numFrames,1); % inside
outsideMasks = cell(numFrames,1); % outside
bwSph = cell(numFrames,1); % mask of the spheroid

for k = 1:numFrames
    
    frame = movieDataBF(:,:,k);
    
    % smooth with Gaussian filter
    frameFilt = imgaussfilt(frame, 5);
    
    % threshold (Otsu threshold)
    level = graythresh(frameFilt);
    bw = frameFilt < level*255;
    
    % clean mask
    bw = bwareaopen(bw, 500); % remove noise
    bw = imfill(bw, 'holes'); % fill inside

    % keep largest object
    bw = bwareafilt(bw, 1);
    bwSph{k} = bw;
    
    % boundary
    B = bwboundaries(bw);
    boundaries{k} = B{1};
    
    % distances
    distOutside = bwdist(bw);
    distInside  = bwdist(~bw);
    
    % around belt (20 outside + 10 inside)
    belt = (distOutside <= 20 & ~bw) | (distInside <= 10 & bw);
    beltMasks{k} = belt;
    
    % inside (shrink object by 10 px)
    inside = distInside > 10;   % keep pixels deeper than 10 px inside
    insideMasks{k} = inside;
    
    % outside (20 px < dist < 60 px)
    outside = distOutside > 20 & distOutside <= 60 & ~bw;
    outsideMasks{k} = outside;
    
end

%% Overlay the boundary, belt, and inside area on the movie

for k = 1:numFrames
    
    imshow(movieDataBF(:,:,k), []);
    hold on;
    
    % around belt
    belt = beltMasks{k};
    [y_belt, x_belt] = find(belt);
    plot(x_belt, y_belt, '.', 'Color', [0.75 0.25 0.75], 'MarkerSize', 0.1); % magenta
    
    % inside
    inside = insideMasks{k};
    [y_in, x_in] = find(inside);
    plot(x_in, y_in, '.', 'Color', [0.00 0.8 0.8], 'MarkerSize', 0.1); % cyan
    
    % outside
    outside = outsideMasks{k};
    [y_out, x_out] = find(outside);
    plot(x_out, y_out, '.', 'Color', [0.0660 0.4430 0.7450], 'MarkerSize', 2); % blue
    
    % boundary
    B = boundaries{k};
    plot(B(:,2), B(:,1), 'Color', [1.0 0.8 0], 'LineWidth', 2); % yellow
    
    hold off;
    pause(0.05);
end

%% Fitting with the circle to find the middle

centers = zeros(numFrames, 2);   % [x, y]
radii   = zeros(numFrames, 1);   % radii
 
for k = 1:numFrames
    
    B = boundaries{k};
    
    x = B(:,2);   % column → x
    y = B(:,1);   % row → y
    
    % circle fit (least squares)
    A = [x, y, ones(size(x))];
    b = -(x.^2 + y.^2);
    
    params = A \ b;
    
    xc = -params(1)/2;
    yc = -params(2)/2;
    r  = sqrt((xc^2 + yc^2) - params(3));
    
    centers(k,:) = [xc, yc];
    radii(k) = r;
end

% visual check
theta = linspace(0, 2*pi, 200);

for k = 1:numFrames
    
    imshow(movieDataBF(:,:,k), []);
    hold on;
    
    % boundary
    B = boundaries{k};
    plot(B(:,2), B(:,1), 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2);
    
    % fitted circle
    xc = centers(k,1);
    yc = centers(k,2);
    r  = radii(k);
    
    xfit = xc + r*cos(theta);
    yfit = yc + r*sin(theta);
    
    plot(xfit, yfit, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 2);
    
    % center point
    plot(xc, yc, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
    
    hold off;
    pause(0.05);
end

%% Fluorescent channels (GFP & mCherry)

infoGFP = imfinfo(filenameGFP);
infomCh = imfinfo(filenamemCh);

movieDataGFP = zeros(infoGFP(1).Height, infoGFP(1).Width, numFrames, 'uint8');
movieDatamCh = zeros(infomCh(1).Height, infomCh(1).Width, numFrames, 'uint8');

for k = 1:numFrames
    movieDataGFP(:,:,k) = imread(filenameGFP, k);
    movieDatamCh(:,:,k) = imread(filenamemCh, k);
end

% Compute total intensity in the belt, inside and outside ring of the spheroid
% for both fluorescent channels

SphArea = zeros(numFrames,1);

% GFP
sumBeltIntensity_GFP   = zeros(numFrames,1);
meanInsideIntensity_GFP = zeros(numFrames,1);
sumInsideIntensity_GFP = zeros(numFrames,1);

meanOutsideIntensity_GFP = zeros(numFrames,1);
sumOutsideIntensity_GFP = zeros(numFrames,1);

% mCherry
sumBeltIntensity_mCh   = zeros(numFrames,1);
meanInsideIntensity_mCh = zeros(numFrames,1);
sumInsideIntensity_mCh = zeros(numFrames,1);

meanOutsideIntensity_mCh = zeros(numFrames,1);
sumOutsideIntensity_mCh = zeros(numFrames,1);

for k = 1:numFrames
    
    % size of the spheroid
    SphArea(k) = sum(inside(:));
    
    % GFP channel
    frameGFP = movieDataGFP(:,:,k);
    
    belt = beltMasks{k};
    inside = insideMasks{k};
    outside = outsideMasks{k};
    
    % belt
    pixelsBelt_GFP = double(frameGFP(belt));
    sumBeltIntensity_GFP(k) = sum(pixelsBelt_GFP);
    
    % inside
    pixelsInside_GFP = double(frameGFP(inside));
    meanInsideIntensity_GFP(k) = mean(pixelsInside_GFP);
    sumInsideIntensity_GFP(k) = sum(pixelsInside_GFP);
    
    % outside
    pixelsOutside_GFP = double(frameGFP(outside));
    meanOutsideIntensity_GFP(k) = mean(pixelsOutside_GFP);
    sumOutsideIntensity_GFP(k) = sum(pixelsOutside_GFP);
    
    % mCherry channel
    frameMCh = movieDatamCh(:,:,k);
    
    % belt
    pixelsBelt_mCh = double(frameMCh(belt));
    sumBeltIntensity_mCh(k) = sum(pixelsBelt_mCh);
    
    % inside
    pixelsInside_mCh = double(frameMCh(inside));
    meanInsideIntensity_mCh(k) = mean(pixelsInside_mCh);
    sumInsideIntensity_mCh(k) = sum(pixelsInside_mCh);
    
    % outside
    pixelsOutside_mCh = double(frameMCh(outside));
    meanOutsideIntensity_mCh(k) = mean(pixelsOutside_mCh);
    sumOutsideIntensity_mCh(k) = sum(pixelsOutside_mCh);
       
end

%% Normalize Intensity (Inside and Around)
% for two fluorescent channels: GFP and mCherry

time = 0:(numFrames-1);
time = timeframeint * time;

% normalize Belt intensity to Belt intensity at t=0
sumBeltIntensity_GFP_norm  = sumBeltIntensity_GFP / sumBeltIntensity_GFP(1);
sumBeltIntensity_mCh_norm  = sumBeltIntensity_mCh / sumBeltIntensity_mCh(1);

% normalize Inside intensity to the average Outside intensity at t=0 by the number
% of pixels of Inside area
sumInsideIntensity_GFP_norm = SphArea.*(meanInsideIntensity_GFP-meanInsideIntensity_GFP(1)) / (meanOutsideIntensity_GFP(1)*SphArea(1));
sumInsideIntensity_GFP_norm(sumInsideIntensity_GFP_norm < 0) = 0; 

sumInsideIntensity_mCh_norm = SphArea.*(meanInsideIntensity_mCh-meanInsideIntensity_mCh(1)) / (meanOutsideIntensity_mCh(1)*SphArea(1));
sumInsideIntensity_mCh_norm(sumInsideIntensity_mCh_norm < 0) = 0; 

% smooth the curves
sumBeltIntensity_GFP_norm_sm  = smoothdata(sumBeltIntensity_GFP_norm, 'movmean', 3);
sumInsideIntensity_GFP_norm_sm = smoothdata(sumInsideIntensity_GFP_norm, 'movmean', 3);
sumBeltIntensity_mCh_norm_sm  = smoothdata(sumBeltIntensity_mCh_norm, 'movmean', 3);
sumInsideIntensity_mCh_norm_sm = smoothdata(sumInsideIntensity_mCh_norm, 'movmean', 3);

%% Plot total GFP Intensity with two y-axes
% both fluorescent channels

figure('Position', [100 100 550 300]);

c1 = [0.90 0.75 0.15]; % GFP channel (yellow)
c2 = [0.75 0.25 0.75]; % mCherry channel (magenta)

% Left axis (around)
% GFP channel
yyaxis left

s1 = scatter(time, sumBeltIntensity_GFP_norm, ...
    35, 'MarkerFaceColor', c1, ...
    'MarkerEdgeColor', c1);
s1.MarkerFaceAlpha = 0.2;
s1.MarkerEdgeAlpha = 0.6;
hold on;

p1 = plot(time, sumBeltIntensity_GFP_norm_sm, ...
    'Color', c1, ...
    'LineWidth', 2, ...
    'LineStyle', '-');

% mCherry channel
s2 = scatter(time, sumBeltIntensity_mCh_norm, ...
    35, 'MarkerFaceColor', c2, ...
    'MarkerEdgeColor', c2);
s2.MarkerFaceAlpha = 0.2;
s2.MarkerEdgeAlpha = 0.6;
hold on;

p2 = plot(time, sumBeltIntensity_mCh_norm_sm, ...
    'Color', c2, ...
    'LineWidth', 2, ...
    'LineStyle', '-');

ylabel('Total Intensity (Around)');

ax = gca;
ax.YAxis(1).Color = [0 0 0];   % left axis color
ylim([0.8 5.5]);

% right axis (Inside)
% GFP channel
yyaxis right

s3 = scatter(time, sumInsideIntensity_GFP_norm, 35, 'd', ...
    'MarkerFaceColor', c1, ...
    'MarkerEdgeColor', c1);
s3.MarkerFaceAlpha = 0.2;
s3.MarkerEdgeAlpha = 0.6;
hold on;

p3 = plot(time, sumInsideIntensity_GFP_norm_sm, ...
    'Color', c1, ...
    'LineWidth', 2, ...
    'LineStyle', '-'); 

% mCherry channel
s4 = scatter(time, sumInsideIntensity_mCh_norm, 35, 'd', ...
    'MarkerFaceColor', c2, ...
    'MarkerEdgeColor', c2);
s4.MarkerFaceAlpha = 0.2;
s4.MarkerEdgeAlpha = 0.6;
hold on;

p4 = plot(time, sumInsideIntensity_mCh_norm_sm, ...
    'Color', c2, ...
    'LineWidth', 2, ...
    'LineStyle', '-');

ylabel('Total Intensity (Inside)');
ax.YAxis(2).Color = [0.45 0.45 0.45];   % right axis color

% align left y=1 with right y=0
yyaxis left
yl = ylim;
frac = (1 - yl(1)) / (yl(2) - yl(1));
yyaxis right
ymaxR = max(sumInsideIntensity_GFP_norm_sm)*1.1;
yminR = -frac*ymaxR/(1-frac);
ylim([yminR ymaxR]);

% reference line
yyaxis left
yline(1,'--','Color',[.35 .35 .35],'LineWidth',1.5);

xlabel('Time, min');
xlim([0 30]);

legend([p1 p2 p3 p4], {'Around GFP', 'Around mCh', 'Inside GFP', 'Inside mCh'});

%% Plot the average (across angles) intensity as a distance from the center
% Radial intensity profile around fitted spheroid center

% maximum radius to analyze (pixels)
maxRadius = 256;

radialProfilesGFP = zeros(numFrames, maxRadius);
radialProfilesmCh = zeros(numFrames, maxRadius);

% image coordinates
[xx, yy] = meshgrid(1:size(movieDataGFP,2), 1:size(movieDataGFP,1));

for k = 1:numFrames
    
    % current fluorescence frame
    frameGFP = double(movieDataGFP(:,:,k));
    framemCh = double(movieDatamCh(:,:,k));
    
    % fitted center
    xc = centers(k,1);
    yc = centers(k,2);
    
    % radial distance map from fitted center
    rr = sqrt((xx - xc).^2 + (yy - yc).^2);
    
    % compute mean intensity in radial bins
    for r = 1:maxRadius
        
        mask = rr >= (r-1) & rr < r;
        pixelsGFP = frameGFP(mask);
        pixelsmCh = framemCh(mask);
        radialProfilesGFP(k,r) = mean(pixelsGFP);
        radialProfilesmCh(k,r) = mean(pixelsmCh);
    end
end

% Smooth radial profiles across radius (5-pixel window)
radialProfilesGFP_sm = zeros(size(radialProfilesGFP));
radialProfilesmCh_sm = zeros(size(radialProfilesmCh));

for k = 1:numFrames
    
    radialProfilesGFP_sm(k,:) = smoothdata(radialProfilesGFP(k,:), ...
        'movmean', 3);
    radialProfilesmCh_sm(k,:) = smoothdata(radialProfilesmCh(k,:), ...
        'movmean', 3);
end

%% Angular scan within Outside mask (average intensity within the sector)
% both fluorescent channels

sectorWidth = 10;
sectorEdges = 0:sectorWidth:360;
sectorCenters = sectorEdges(1:end-1) + sectorWidth/2;

numSectors = length(sectorCenters);

angleTimeMapGFP = nan(numFrames, numSectors);
angleTimeMapmCh = nan(numFrames, numSectors);

[xx,yy] = meshgrid(1:size(movieDataGFP,2), ...
                   1:size(movieDataGFP,1));

for k = 1:numFrames

    frameGFP = double(movieDataGFP(:,:,k));
    framemCh = double(movieDatamCh(:,:,k));

    xc = centers(k,1);
    yc = centers(k,2);

    outside = outsideMasks{k};

    % polar coordinates
    theta = atan2d(xx-xc, -(yy-yc));   % 0° = top
    theta(theta<0) = theta(theta<0)+360;

    for s = 1:numSectors

        alpha1 = sectorEdges(s);
        alpha2 = sectorEdges(s+1);

        sectorMask = theta >= alpha1 & theta < alpha2;

        mask = outside & sectorMask;

        pixelsGFP = frameGFP(mask);
        pixelsmCh = framemCh(mask);

        if ~isempty(pixelsGFP)
            angleTimeMapGFP(k,s) = mean(pixelsGFP);
        end
        
        if ~isempty(pixelsmCh)
            angleTimeMapmCh(k,s) = mean(pixelsmCh);
        end

    end
end

%% Plot average intensity within angular sectors VS Time
figure('Position',[100 100 565 700])

% GFP
ax1 = subplot(2,1,1);
hold(ax1,'on')

imagesc(time, sectorCenters, angleTimeMapGFP')

hold(ax1,'off')

axis xy

xlim([0 30]);
ylim([0 360]);
xlabel('Time [min]')
ylabel('Angle [deg]')

colormap(ax1, parula)

cb1 = colorbar(ax1);
cb1.Label.String = 'Average GFP Intensity';

% mCh
ax2 = subplot(2,1,2);
hold(ax2,'on')

imagesc(time, sectorCenters, angleTimeMapmCh')

hold(ax2,'off')

axis xy

xlim([0 30]);
ylim([0 360]);
xlabel('Time [min]')
ylabel('Angle [deg]')

colormap(ax2, cool)

cb2 = colorbar(ax2);
cb2.Label.String = 'Average mCherry Intensity';

box(ax2,'on')


%% Angular scan within Outside mask as a function of distance from spheroid
% both fluorescent channels

sectorWidth = 10;
sectorEdges = 0:sectorWidth:360;
sectorCenters = sectorEdges(1:end-1) + sectorWidth/2;
numSectors = length(sectorCenters);

distanceBins = 20:1:60;      % px from spheroid boundary
numDistBins = length(distanceBins)-1;

distanceAngleTimeMapGFP = nan(numFrames, numDistBins, numSectors);
distanceAngleTimeMapmCh = nan(numFrames, numDistBins, numSectors);

[xx,yy] = meshgrid(1:size(movieDataGFP,2), ...
                   1:size(movieDataGFP,1));

for k = 1:numFrames

    frameGFP = double(movieDataGFP(:,:,k));
    framemCh = double(movieDatamCh(:,:,k));

    xc = centers(k,1);
    yc = centers(k,2);

    outside = outsideMasks{k};
    spheroidMask = bwSph{k};
    
    % distance from spheroid boundary
    distMap = bwdist(spheroidMask);

    % polar coordinates
    theta = atan2d(xx-xc, -(yy-yc));   % 0° = top
    theta(theta<0) = theta(theta<0)+360;

    for s = 1:numSectors

        alpha1 = sectorEdges(s);
        alpha2 = sectorEdges(s+1);

        sectorMask = theta >= alpha1 & theta < alpha2;

        for d = 1:numDistBins

            d0 = distanceBins(d);
            d1 = distanceBins(d+1);

            distanceMask = distMap >= d0 & distMap < d1;

            mask = outside & sectorMask & distanceMask;

            pixelsGFP = frameGFP(mask);
            pixelsmCh = framemCh(mask);

            if ~isempty(pixelsGFP)
                distanceAngleTimeMapGFP(k,d,s) = mean(pixelsGFP);
            end
            
            if ~isempty(pixelsmCh)
                distanceAngleTimeMapmCh(k,d,s) = mean(pixelsmCh);
            end
        end
    end
end

%% Distance-time map for a selected angular region
% expect manually where the wave is stronger and synchronized
% for both fluorescent channels

alpha1 = 0;   % degrees; MANUALLY SELECT THE ANGLE
alpha2 = 360;   % degrees; MANUALLY SELECT THE ANGLE

distanceBins = 20:1:60;      % px from spheroid boundary
numDistBins = length(distanceBins)-1;

% time × distance
distanceTimeMapGFP = nan(numFrames, numDistBins);
distanceTimeMapmCh = nan(numFrames, numDistBins);

[xx,yy] = meshgrid(1:size(movieDataGFP,2), ...
                   1:size(movieDataGFP,1));

for k = 1:numFrames

    frameGFP = double(movieDataGFP(:,:,k));
    framemCh = double(movieDatamCh(:,:,k));

    xc = centers(k,1);
    yc = centers(k,2);

    outside = outsideMasks{k};
    spheroidMask = bwSph{k};

    % distance from spheroid boundary
    distMap = bwdist(spheroidMask);

    % polar coordinates
    theta = atan2d(xx-xc, -(yy-yc));   % 0° = top
    theta(theta<0) = theta(theta<0)+360;

    % angular mask
    if alpha1 < alpha2
        sectorMask = theta >= alpha1 & theta < alpha2;
    else
        % crossing 0°
        sectorMask = theta >= alpha1 | theta < alpha2;
    end

    for d = 1:numDistBins

        d0 = distanceBins(d);
        d1 = distanceBins(d+1);

        distanceMask = distMap >= d0 & distMap < d1;

        mask = outside & sectorMask & distanceMask;

        pixelsGFP = frameGFP(mask);
        pixelsmCh = framemCh(mask);

        if ~isempty(pixelsGFP)
            distanceTimeMapGFP(k,d) = mean(pixelsGFP);
        end
        
        if ~isempty(pixelsmCh)
            distanceTimeMapmCh(k,d) = mean(pixelsmCh);
        end
        
    end
end

%% Plot distance-time map for a selected angular region
% for both channels

figure('Position',[100 100 300 600])

% GFP
ax1 = subplot(2,1,1);
hold(ax1,'on')

imagesc(time,...
        pixelsize*distanceBins(1:end-1),...
        distanceTimeMapGFP')

axis xy

ylim([54 140]);
xlabel('Time [min]')
ylabel('Distance from spheroid edge [\mum]');

hold(ax1,'off')

colormap(ax1, parula)
cb1 = colorbar(ax1);
cb1.Label.String = 'Average GFP Intensity';

box(ax1,'on')

% mCh
ax2 = subplot(2,1,2);
hold(ax2,'on')

imagesc(time,...
        pixelsize*distanceBins(1:end-1),...
        distanceTimeMapmCh')

hold(ax2,'off')
axis xy

ylim([54 140]);
xlabel('Time [min]')
ylabel('Distance from spheroid edge [\mum]');

hold(ax1,'off')

colormap(ax2, cool)
cb2 = colorbar(ax2);
cb2.Label.String = 'Average mCherry Intensity';

box(ax2,'on')

%% Total intensity in a selected angular range
% both fluorescent channels

sectorIntensityGFP = nan(numFrames,1);
sectorIntensitymCh = nan(numFrames,1);

[xx,yy] = meshgrid(1:size(movieDataGFP,2), ...
                   1:size(movieDataGFP,1));

for k = 1:numFrames

    frameGFP = double(movieDataGFP(:,:,k));
    framemCh = double(movieDatamCh(:,:,k));

    xc = centers(k,1);
    yc = centers(k,2);

    % angle map
    theta = atan2d(xx-xc, -(yy-yc));   % 0° = top
    theta(theta<0) = theta(theta<0)+360;

    % angular mask
    if alpha1 < alpha2
        sectorMask = theta >= alpha1 & theta < alpha2;
    else
        % crossing 0°
        sectorMask = theta >= alpha1 | theta < alpha2;
    end

    % outside region
    outside = outsideMasks{k};
    spheroidMask = bwSph{k};

    % distance from spheroid boundary
    distMap = bwdist(spheroidMask);

    d0 = 20; % px
    d1 = 60; % px
    distanceMask = distMap >= d0 & distMap < d1;
    
    % combined mask
    mask = outside & sectorMask & distanceMask;

    pixelsGFP = frameGFP(mask);
    if ~isempty(pixelsGFP)
        sectorIntensityGFP(k) = sum(pixelsGFP);
    end

    pixelsmCh = framemCh(mask);
    if ~isempty(pixelsmCh)
        sectorIntensitymCh(k) = sum(pixelsmCh);
    end
    
end

% normalize total intensity to t=0
sectorIntensityGFP_norm = sectorIntensityGFP / sectorIntensityGFP(1);
sectorIntensitymCh_norm = sectorIntensitymCh / sectorIntensitymCh(1);

sectorIntensityGFP_norm_sm = smoothdata( ...
    sectorIntensityGFP_norm, ...
    'movmean', ...
    3);

sectorIntensitymCh_norm_sm = smoothdata( ...
    sectorIntensitymCh_norm, ...
    'movmean', ...
    3);

%% Plot Total intensity (outside ring, selected angle range)

figure('Position',[100 100 550 300])

plot(time, sectorIntensityGFP_norm_sm, '-', ...
    'LineWidth',2, ...
    'Color',c1); hold on;

plot(time, sectorIntensitymCh_norm_sm, '-', ...
    'LineWidth',2, ...
    'Color',c2);

xlabel('Time [min]')
ylabel('Total Intensity [AU]');
title(sprintf('Sector %d°-%d°', alpha1, alpha2))
box off;


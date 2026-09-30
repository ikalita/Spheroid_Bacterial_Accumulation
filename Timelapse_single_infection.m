%% Image analysis of bacterial accumulation around spheroids from timelapse sequence
%
% This script quantifies the spatial and temporal dynamics of bacterial
% accumulation around spheroids from time-lapse microscopy images.
%
% Workflow:
% 1) Identify the spheroid boundary in the bright-field channel.
% 2) Quantify fluorescence intensity around and inside the spheroid.
% 3) Generate radial fluorescence intensity profiles over time.
% 4) Identify angular sectors with pronounced dispersal wave
%
% The spheroid position is identified independently in each frame,
% and the spheroid center is estimated by fitting a circle to its boundary.

clc
clear all
close all

% microscopy details
pixelsize = 2.77; % in microns
timeframeint = 30/60; % time between frames in minutes

c1 = [0 0 0]; % Surface (black)
c2 = [0.5 0.5 0.5]; % Inside (grey)

% load the bright-field movie
filenameBF = 'example_BF.tif'; % CHANGE FILE NAME (bright-field channel)
% load the fluorescent movie
filenameGFP = 'example_GFP.tif'; % CHANGE FILE NAME (fluorescent channel)

infoBF = imfinfo(filenameBF);
numFrames = numel(infoBF);

movieDataBF = zeros(infoBF(1).Height, infoBF(1).Width, numFrames, 'uint8');

for k = 1:numFrames
    movieDataBF(:,:,k) = imread(filenameBF, k);
end

%% Identify the position of the spheroid for each frame

boundaries = cell(numFrames,1);
beltMasks = cell(numFrames,1); % surface
insideMasks = cell(numFrames,1); % inside
outsideMasks = cell(numFrames,1); % all outside rings 
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
    
    % surface (20 outside + 10 inside)
    belt = (distOutside <= 20 & ~bw) | (distInside <= 10 & bw);
    beltMasks{k} = belt;
    
    % inside (shrink object by 10 px)
    inside = distInside > 10;   % keep pixels deeper than 10 px inside
    insideMasks{k} = inside;
    
    % All rings together: (20 px < dist <= 60 px)
    outside = distOutside > 20 & distOutside <= 60 & ~bw;
    outsideMasks{k} = outside;
    
end

%% Overlay the boundary, surface, and inside area on the movie

for k = 1:numFrames
    
    imshow(movieDataBF(:,:,k), []);
    hold on;
    
    % boundary
    B = boundaries{k};
    plot(B(:,2), B(:,1), 'Color', [0, 0.4470, 0.7410], 'LineWidth', 4); % blue
    
    % surface
    belt = beltMasks{k};
    [y_belt, x_belt] = find(belt);
    plot(x_belt, y_belt, '.', 'Color', [0.8500, 0.3250, 0.0980], 'MarkerSize', 2); % orange
    
    % inside
    inside = insideMasks{k};
    [y_in, x_in] = find(inside);
    plot(x_in, y_in, '.', 'Color', [0.4660, 0.6740, 0.1880], 'MarkerSize', 2); % green
    
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

%% Fluorescent channel

infoGFP = imfinfo(filenameGFP);
movieDataGFP = zeros(infoGFP(1).Height, infoGFP(1).Width, numFrames, 'uint8');

for k = 1:numFrames
    movieDataGFP(:,:,k) = imread(filenameGFP, k);
end

% Compute total intensity in the surface, inside and outside rings of the spheroid

sumBeltIntensity   = zeros(numFrames,1);
meanInsideIntensity = zeros(numFrames,1);
sumInsideIntensity = zeros(numFrames,1);
SphArea = zeros(numFrames,1);
meanOutsideIntensity = zeros(numFrames,1);

for k = 1:numFrames
    
    frameFM = movieDataGFP(:,:,k);
    
    % belt
    belt = beltMasks{k};
    pixelsBelt = double(frameFM(belt));
    sumBeltIntensity(k) = sum(pixelsBelt);
    
    % inside
    inside = insideMasks{k};
    pixelsInside = double(frameFM(inside));
    meanInsideIntensity(k) = mean(pixelsInside);
    sumInsideIntensity(k) = sum(pixelsInside);
    SphArea(k) = sum(inside(:));
    
    % outside
    outside = outsideMasks{k};
    pixelsOutside = double(frameFM(outside));
    meanOutsideIntensity(k) = mean(pixelsOutside);
    
end

%% Normalize Intensity (Inside and Around)

time = 0:(numFrames-1);
time = timeframeint * time;

% normalize Belt intensity to Belt intensity at t=0
sumBeltIntensity_norm  = sumBeltIntensity / sumBeltIntensity(1);

% normalize Inside intensity to the average Outside intensity at t=0 by the number
% of pixels of Inside area
sumInsideIntensity_norm = SphArea.*(meanInsideIntensity-meanInsideIntensity(1)) / (meanOutsideIntensity(1)*SphArea(1));
sumInsideIntensity_norm(sumInsideIntensity_norm < 0) = 0; 

% smooth the curves
sumBeltIntensity_norm_sm  = smoothdata(sumBeltIntensity_norm, 'movmean', 3);
sumInsideIntensity_norm_sm = smoothdata(sumInsideIntensity_norm, 'movmean', 3);


%% Plot total GFP Intensity with two y-axes

figure('Position', [100 100 550 300]);

% Left axis (around)
yyaxis left

s1 = scatter(time, sumBeltIntensity_norm, ...
    50, 'MarkerFaceColor', c1, ...
    'MarkerEdgeColor', c1);

s1.MarkerFaceAlpha = 0.45;
s1.MarkerEdgeAlpha = 0.85;
hold on;

p1 = plot(time, sumBeltIntensity_norm_sm, ...
    'Color', c1, ...
    'LineWidth', 2, ...
    'LineStyle', '-');
ylabel('Total GFP Intensity (Around)');

ax = gca;
ax.YAxis(1).Color = c1;   % left axis color
ylim([0.8 3.5]);

% right axis (Inside)
yyaxis right

s2 = scatter(time, sumInsideIntensity_norm, ...
    50, 'd', ...
    'MarkerFaceColor', c2, ...
    'MarkerEdgeColor', c2);

s2.MarkerFaceAlpha = 0.15;
s2.MarkerEdgeAlpha = 0.5;
hold on;

p2 = plot(time, sumInsideIntensity_norm_sm, ...
    'Color', c2, ...
    'LineWidth', 2, ...
    'LineStyle', '-'); 

ylabel('Total GFP Intensity (Inside)');
ax.YAxis(2).Color = c2;   % right axis color

% align left y=1 with right y=0
yyaxis left
yl = ylim;
frac = (1 - yl(1)) / (yl(2) - yl(1));
yyaxis right
ymaxR = max(sumInsideIntensity_norm_sm)*1.1;
yminR = -frac*ymaxR/(1-frac);
ylim([yminR ymaxR]);

% reference line
yyaxis left
yline(1,'--','Color',[.5 .5 .5],'LineWidth',1.5);

xlabel('Time, min');
xlim([0 45]);

legend([p1 p2], {'Surface', 'Interior'});

box on;

%% Plot the average (across angles) intensity as a distance from the center
% Radial intensity profile around fitted spheroid center

% maximum radius to analyze (pixels)
maxRadius = 256;

radialProfiles = zeros(numFrames, maxRadius);

% image coordinates
[xx, yy] = meshgrid(1:size(movieDataGFP,2), 1:size(movieDataGFP,1));

for k = 1:numFrames
    
    % current fluorescence frame
    frameFM = double(movieDataGFP(:,:,k));
    
    % fitted center
    xc = centers(k,1);
    yc = centers(k,2);
    
    % radial distance map from fitted center
    rr = sqrt((xx - xc).^2 + (yy - yc).^2);
    
    % compute mean intensity in radial bins
    for r = 1:maxRadius
        
        mask = rr >= (r-1) & rr < r;
        pixels = frameFM(mask);
        radialProfiles(k,r) = mean(pixels);
    end
end

% Smooth radial profiles across radius (5-pixel window)
radialProfiles_sm = zeros(size(radialProfiles));

for k = 1:numFrames
    
    radialProfiles_sm(k,:) = smoothdata(radialProfiles(k,:), ...
        'movmean', 5);
end

%% Angular scan within Outside mask (average intensity within the sector)

sectorWidth = 10;
sectorEdges = 0:sectorWidth:360;
sectorCenters = sectorEdges(1:end-1) + sectorWidth/2;

numSectors = length(sectorCenters);

angleTimeMap = nan(numFrames, numSectors);

[xx,yy] = meshgrid(1:size(movieDataGFP,2), ...
                   1:size(movieDataGFP,1));

for k = 1:numFrames

    frameFM = double(movieDataGFP(:,:,k));

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

        pixels = frameFM(mask);

        if ~isempty(pixels)
            angleTimeMap(k,s) = mean(pixels);
        end

    end
end

%% Plot average intensity within angular sectors VS Time
figure('Position',[100 100 550 300])

imagesc(time, sectorCenters, angleTimeMap')

axis xy

xlim([0 15]);
xlabel('Time [min]')
ylabel('Angle [deg]')

title('Intensity within angular sectors VS Time')
colorbar

%% Distance-time map for a selected angular region
% expect manually where the wave is stronger and synchronized

alpha1 = 0;   % degrees; MANUALLY SELECT THE ANGLE
alpha2 = 360;   % degrees; MANUALLY SELECT THE ANGLE

distanceBins = 20:1:60;      % px from spheroid boundary
numDistBins = length(distanceBins)-1;

% time × distance
distanceTimeMap = nan(numFrames, numDistBins);

[xx,yy] = meshgrid(1:size(movieDataGFP,2), ...
                   1:size(movieDataGFP,1));

for k = 1:numFrames

    frameFM = double(movieDataGFP(:,:,k));

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

        pixels = frameFM(mask);

        if ~isempty(pixels)
            distanceTimeMap(k,d) = mean(pixels);
        end
    end
end

%% Plot distance-time map for a selected angular region

figure('Position',[100 100 550 300])

imagesc(time,...
        pixelsize*distanceBins(1:end-1),...
        distanceTimeMap')

axis xy

xlabel('Time [min]')
ylabel('Distance from spheroid edge [\mum]');

title(sprintf('%d°-%d°', alpha1, alpha2))

colorbar

xlim([0 20]);
ylim([54 140]);

%% Average intensity (outside rings) vs Time in a selected angular range

sectorIntensity = nan(numFrames,1);

[xx,yy] = meshgrid(1:size(movieDataGFP,2), ...
                   1:size(movieDataGFP,1));

for k = 1:numFrames

    frameFM = double(movieDataGFP(:,:,k));

    xc = centers(k,1);
    yc = centers(k,2);

    % angle map
    theta = atan2d(xx-xc, -(yy-yc));   % 0° = top
    theta(theta<0) = theta(theta<0)+360;

    % selected angular sector
    sectorMask = theta >= alpha1 & theta < alpha2;

    % outside rings
    outside = outsideMasks{k};
    
    % combined mask
    mask = outside & sectorMask;
    pixels = frameFM(mask);
    
    if ~isempty(pixels)
        sectorIntensity(k) = sum(pixels);
    end

end

% normalize total intensity to t=0
sectorIntensity_norm = sectorIntensity / sectorIntensity(1);

sectorIntensity_norm_sm = smoothdata( ...
    sectorIntensity_norm, ...
    'movmean', ...
    3);

%% Plot Total intensity (outside ring, selected angle range)

figure('Position',[100 100 550 300])

plot(time,sectorIntensity_norm_sm,'LineWidth',2); 

ylabel('Total GFP intensity [AU]')
xlabel('Time [min]')

title(sprintf('Sector %d°-%d°', alpha1, alpha2))

box on

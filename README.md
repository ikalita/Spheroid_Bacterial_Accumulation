# Bacterial Accumulation Dynamics Around Spheroids

MATLAB script for analyzing the spatial and temporal dynamics of bacterial fluorescence around spheroids from time-lapse microscopy images.

## Workflow

The analysis consists of the following steps:

1. Identify the spheroid boundary in the bright-field channel.
2. Define surface, inside, and outside regions relative to the spheroid.
3. Quantify fluorescence intensity inside and around the spheroid over time.
4. Generate radial fluorescence intensity profiles to characterize bacterial accumulation as a function of distance from the spheroid.
5. Generate angular intensity maps to identify sectors with pronounced fluorescence accumulation.
6. Generate distance–time maps for a selected angular sector to visualize the spatial and temporal dynamics of bacterial accumulation.

## Input

The script requires:

* Bright-field time-lapse TIFF movie
* Fluorescence time-lapse TIFF movie (e.g., GFP)
* Pixel size
* Frame interval

Update the input filenames and analysis parameters at the beginning of the script before running the analysis.

## Output

The script generates:

* Fluorescence intensity plots over time
* Radial fluorescence intensity profiles
* Angular intensity maps
* Distance–time maps for selected angular regions

## Requirements

* MATLAB
* Image Processing Toolbox
* Bright-field time-lapse TIFF images
* Fluorescence time-lapse TIFF images

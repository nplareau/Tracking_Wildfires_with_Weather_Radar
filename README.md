# Tracking_Wildfires_with_Weather_Radar

This directory contains codes that support "Tracking Wildfires with Weather Radars"

Lareau, N. P., Donohoe, A., Roberts, M., & Ebrahimian H., (2022): Tracking Wildfires with Weather Radars. Journal of Geophysical Research - Atmospheres. 

The directory contains two MATLAB codes
(1) NEXRAD_PREPROCESSOR.m: This code reads netcdf NEXRAD radar files, which are obtained using the NOAA Weather and Climate Toolkit (https://www.ncdc.noaa.gov/wct/)
The code outputs a ".mat" file with the relevant radar data and spatial grids. These grids are then used in the next script. 

(2) radar_perimeters_v1.m: This code estimates fire perimeters from the preprocesses NEXRAD data. Details of the code and approach are found in the journal article and in the header of the code

questions about the code should be sent to Neil P. Lareau at nlareau@unr.edu


https://zenodo.org/badge/latestdoi/483526154

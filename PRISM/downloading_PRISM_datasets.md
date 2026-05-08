# Downloading climate data for the on-field experiment. 

Downloading data from PRISM for the US. Varibles from "an" (all networks, refers to time series focused on
providing the best estimate possible, at the expense of temporal consistency) include : precipitation (ppt), minimum temperature (tmin),
maximum temperature (tmax), mean dew point (tdmean), minimum vapor pressure deficit (vpdmin), maximum
vapor pressure deficit (vpdmax), total global shortwave solar irradiance on a horizontal surface (soltotal), and
total global solar irradiance on a sloped surface (solslope). Mean temperature (tmean) is derived as the average of
tmax and tmin. 


Other dataset we could use if not satisfied with PRISM: Daymet https://prism.oregonstate.edu/explorer/ (was >3TB to download daily reports for small spatial scale resolution because I could not select a specific year and it goes back to 1950). 


Using Cyberduck to download daily climate data with small (800m) spatial scale resolution: https://prism.oregonstate.edu/documents/PRISM_downloads_FTP.pdf


Details of the stations, models, variable descriptions: https://prism.oregonstate.edu/documents/PRISM_datasets.pdf 

Once the data is downloaded, unzip the folders to access the raster file (.tif), and remove the zipped folder using:
```
for f in *.zip; do unzip "$f"; rm "$f"; done
```

To extract the data for each day of 2024, from April to October (included) and create a separate dataset for each variable, I used a custom R script that leverages the package "terra" for cog files: 
[extract_climate_data_prism.R](https://github.com/rozenn-pineau/on-field/tree/main/PRISM#:~:text=yesterday-,extract_climate_data_prism.R,-Create%20extract_climate_data_prism.R)


I used a second script to extract climate data relevant to the "in betweens" of our sampling points. For now, this script brings together the climate variables and extract the min, mean and max temperature of each sanpling interval








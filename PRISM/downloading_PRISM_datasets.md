# Downloading climate data for the on-field experiment. 

Other dataset we could use if not satisfied with PRISM: Daymet https://prism.oregonstate.edu/explorer/ (was >3TB to download daily reports for small spatial scale resolution because I could not select a specific year and it goes back to 1950). 


Using Cyberduck to download daily climate data with small (800m) spatial scale resolution: https://prism.oregonstate.edu/documents/PRISM_downloads_FTP.pdf


Details of the stations, models, variable descriptions: https://prism.oregonstate.edu/documents/PRISM_datasets.pdf 

Once the data is downloaded, unzip the folders to access the raster file (.tif), and remove the zipped folder using:
```
for f in *.zip; do unzip "$f"; rm "$f"; done
```

To extract the data for each day of 2024 and create a separate dataset for each variable, I used a custom R script that leverages the package "terra" for cog files (example here for precipitation):

```
#R

# libraries
library(colourvalues)
library(extrafont)
library(terra) #cog package

#sampling location coordinates
sampling <- fread("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/metadata/Sample_info.csv")
head(sampling)
#keep one line per sampling site (keeping climate data for the year)
sampling_unique <- sampling %>% group_by(PopEnv) %>% slice_head()
#define locations to extract climate data: convert to SpatVector (points) - create a pointer
points <- vect(sampling_unique, geom=c("long", "lat"), crs="EPSG:4326")
#The primary EPSG code for the NAD83 (North American Datum of 1983) geographic coordinate system is EPSG:4269
n_sampling_points <- dim(points)[1]

# ----------------------------------------------------------------------------- #
# PRCP
# ----------------------------------------------------------------------------- #
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/ppt/")
files <- list.files(pattern = "\\.tif$", full.names = TRUE)
n_days <- length(files)
prcp <- matrix(NA, n_sampling_points, n_days)

for (day in 1:n_days){
  
  # Open the COG
  f <- files[day]
  r <- rast(f)
  print(f)
  
  # Get coordinates (x, y) for all cells
  coords <- crds(r, df=TRUE)
  
  # Extract data
  climate_values <- extract(r, points, ID = F)
  
  # Add to dataset
  prcp[,day] <- climate_values[,1]

}

#prep data for saving
days <- matrix(NA, 1, n_days)
for (i in 1:n_days) {
  days[i] <- strsplit(strsplit(files[i], "_")[[1]][5], ".tif")[[1]][1]}
colnames(prcp) <- days
write.table(prcp, "/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/extracted/ppt_2024.csv", 
            sep = ",", col.names = T, row.names = F, quote = F)

```

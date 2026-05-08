#Goal: I downloaded the daily reports for a selection of climate variables for the 2024 year, for all of North America. 
#These files are in COG format: I need to extract the information from them. 
#I also need to only extract the location, or radius that we are interested in, we don't need the whole continent. 
rm(list= ls())

# libraries
library(colourvalues)
library(extrafont)
library(terra) #cog package
library(data.table)
library(stringr)
library(lubridate)
library(dplyr)

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
# TMEAN
# ----------------------------------------------------------------------------- #

setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/tdmean/")
#keep the weather data for the months of April to October (included)
files <- list.files(pattern = "\\.tif$", full.names = TRUE)
files <- files[month(ymd(str_extract(files, "\\d{8}"))) %in% 4:10]
n_days <- length(files)
tmean <- matrix(NA, n_sampling_points, n_days)

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
  tmean[,day] <- climate_values[,1]
  
}

#Add days as header
days <- matrix(NA, 1, n_days)
for (i in 1:n_days) {
  days[i] <- strsplit(strsplit(files[i], "_")[[1]][5], ".tif")[[1]][1]}
colnames(tmean) <- days
#Add a column to keep location name and coordinates
to_export <- bind_cols(sampling_unique[,c(1:4,8:9)], tmean)
#export data
write.table(to_export, "/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/extracted/tmean_April_October_2024.csv", 
            sep = ",", col.names = T, row.names = F, quote = F)


# ----------------------------------------------------------------------------- #
# TMAX
# ----------------------------------------------------------------------------- #

setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/tmax/")
#keep the weather data for the months of April to October (included)
files <- list.files(pattern = "\\.tif$", full.names = TRUE)
files <- files[month(ymd(str_extract(files, "\\d{8}"))) %in% 4:10]
n_days <- length(files)
tmax <- matrix(NA, n_sampling_points, n_days)

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
  tmax[,day] <- climate_values[,1]
  
}

#Add days as header
days <- matrix(NA, 1, n_days)
for (i in 1:n_days) {
  days[i] <- strsplit(strsplit(files[i], "_")[[1]][5], ".tif")[[1]][1]}
colnames(tmax) <- days
#Add a column to keep location name and coordinates
to_export <- bind_cols(sampling_unique[,c(1:4,8:9)], tmax)
#export data
write.table(to_export, "/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/extracted/tmax_April_October_2024.csv", 
            sep = ",", col.names = T, row.names = F, quote = F)


# ----------------------------------------------------------------------------- #
# TMIN
# ----------------------------------------------------------------------------- #

setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/tmin/")
#keep the weather data for the months of April to October (included)
files <- list.files(pattern = "\\.tif$", full.names = TRUE)
files <- files[month(ymd(str_extract(files, "\\d{8}"))) %in% 4:10]
n_days <- length(files)
tmin <- matrix(NA, n_sampling_points, n_days)

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
  tmin[,day] <- climate_values[,1]
  
}

#Add days as header
days <- matrix(NA, 1, n_days)
for (i in 1:n_days) {
  days[i] <- strsplit(strsplit(files[i], "_")[[1]][5], ".tif")[[1]][1]}
colnames(tmin) <- days
#Add a column to keep location name and coordinates
to_export <- bind_cols(sampling_unique[,c(1:4,8:9)], tmin)
#export data
write.table(to_export, "/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/extracted/tmin_April_October_2024.csv", 
            sep = ",", col.names = T, row.names = F, quote = F)




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

#Add days as header
days <- matrix(NA, 1, n_days)
for (i in 1:n_days) {
  days[i] <- strsplit(strsplit(files[i], "_")[[1]][5], ".tif")[[1]][1]}
colnames(prcp) <- days
#Add a column to keep location name and coordinates
to_export <- bind_cols(sampling_unique[,c(1:4,8:9)], prcp)
#export data
write.table(to_export, "/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/extracted/ppt_2024.csv", 
            sep = ",", col.names = T, row.names = F, quote = F)



# ----------------------------------------------------------------------------- #
# SOLTOT
# ----------------------------------------------------------------------------- #
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/soltot/")
files <- list.files(pattern = "\\.tif$", full.names = TRUE)
files <- files[month(ymd(str_extract(files, "\\d{8}"))) %in% 4:10]
n_days <- length(files)
sol <- matrix(NA, n_sampling_points, n_days)

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
  sol[,day] <- climate_values[,1]
  
}

#Add days as header
days <- matrix(NA, 1, n_days)
for (i in 1:n_days) {
  days[i] <- strsplit(strsplit(files[i], "_")[[1]][5], ".tif")[[1]][1]}
colnames(sol) <- days
#Add a column to keep location name and coordinates
to_export <- bind_cols(sampling_unique[,c(1:4,8:9)], sol)
#export data
write.table(to_export, "/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/extracted/sol_April_October_2024.csv", 
            sep = ",", col.names = T, row.names = F, quote = F)

# ----------------------------------------------------------------------------- #
# VPDMAX
# ----------------------------------------------------------------------------- #
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/vpdmax/")
files <- list.files(pattern = "\\.tif$", full.names = TRUE)
files <- files[month(ymd(str_extract(files, "\\d{8}"))) %in% 4:10]
n_days <- length(files)
vpdmax <- matrix(NA, n_sampling_points, n_days)

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
  vpdmax[,day] <- climate_values[,1]
  
}

#Add days as header
days <- matrix(NA, 1, n_days)
for (i in 1:n_days) {
  days[i] <- strsplit(strsplit(files[i], "_")[[1]][5], ".tif")[[1]][1]}
colnames(vpdmax) <- days
#Add a column to keep location name and coordinates
to_export <- bind_cols(sampling_unique[,c(1:4,8:9)], vpdmax)
#export data
write.table(to_export, "/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/extracted/vpdmax_April_October_2024.csv", 
            sep = ",", col.names = T, row.names = F, quote = F)



#Goal: to extract and summarize the climate data around sampling times.
#We have climate information for every day during that year. 
#We sampled once every month in the summer. We are interested in the climate during the "time in between". 
rm(list= ls())
# libraries
library(colourvalues)
library(extrafont)
library(data.table)
library(stringr)
library(lubridate) #working with dates
library(dplyr)
library(tidyr)
library(purrr)

# -----------------------------
# 1. Load metadata
# -----------------------------
meta <- fread("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/metadata/Sample_info.csv")

meta <- meta %>% mutate(Sample_Date = mdy(Sample_Collection_Date),
  # create matching date format
  Date_ID = format(Sample_Date, "%Y%m%d"))

# -----------------------------
# 2. Load climate datasets
# -----------------------------
ppt <- fread("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/extracted/ppt_April_October_2024.csv")
tmin <- fread("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/extracted/tmin_April_October_2024.csv")
tmean <- fread("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/extracted/tmean_April_October_2024.csv")
tmax <- fread("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/extracted/tmax_April_October_2024.csv")
sol <- fread("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/extracted/sol_April_October_2024.csv")
vpdmax <- fread("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/extracted/vpdmax_April_October_2024.csv")


# -----------------------------
# 3. Convert wide daily columns -> long format
# -----------------------------
ppt_long <- ppt %>%
   pivot_longer(
    cols = matches("^2024"),
    names_to = "Date_ID",
    values_to = "ppt")

tmin_long <- tmin %>%
  pivot_longer(
    cols = matches("^2024"),
    names_to = "Date_ID",
    values_to = "tmin")

tmean_long <- tmean %>%
  pivot_longer(
    cols = matches("^2024"),
    names_to = "Date_ID",
    values_to = "tmean") 

tmax_long <- tmax %>%
  pivot_longer(
    cols = matches("^2024"),
    names_to = "Date_ID",
    values_to = "tmax")

sol_long <- sol %>%
  pivot_longer(
    cols = matches("^2024"),
    names_to = "Date_ID",
    values_to = "sol")

vpdmax_long <- vpdmax %>%
  pivot_longer(
    cols = matches("^2024"),
    names_to = "Date_ID",
    values_to = "vpdmax") 

# -----------------------------
# 4. Join climate to metadata
# -----------------------------

meta_climate <- meta %>% left_join(ppt_long, by = c("PopEnv", "Date_ID"))

head(meta_climate)

climate_summary <- meta %>%
  left_join(ppt_long %>% select(c("PopEnv", "Date_ID"), ppt), by = c("PopEnv", "Date_ID")) %>%
  left_join(tmin_long %>% select(c("PopEnv", "Date_ID"), tmin), by = c("PopEnv", "Date_ID")) %>%
  left_join(tmean_long %>% select(c("PopEnv", "Date_ID"), tmean), by = c("PopEnv", "Date_ID")) %>%
  left_join(tmax_long %>% select(c("PopEnv", "Date_ID"), tmax), by = c("PopEnv", "Date_ID")) %>%
  left_join(sol_long %>% select(c("PopEnv", "Date_ID"), sol), by = c("PopEnv", "Date_ID")) %>%
  left_join(vpdmax_long %>% select(c("PopEnv", "Date_ID"), vpdmax), by = c("PopEnv", "Date_ID")) 

# -----------------------------
# 5. Check
# -----------------------------
#Run a few checks: same location through time
plot(climate_summary$Date_ID, climate_summary$tmean, pch = 16)
plot(climate_summary$Date_ID, climate_summary$sol, pch = 16)
plot(climate_summary$Date_ID, climate_summary$ppt, pch = 16)

#for a single location 
idx <- climate_summary$PopEnv == "Monee1ag"
plot(climate_summary$Date_ID[idx], climate_summary$tmean[idx], pch = 16)
plot(climate_summary$Date_ID[idx], climate_summary$sol[idx], pch = 16)
plot(climate_summary$Date_ID[idx], climate_summary$ppt[idx], pch = 16)

# -----------------------------
# 6. Export
# -----------------------------
write.table(climate_summary, 
            "/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/summaries/climate_time_points.csv", 
            sep = ",", col.names = T, row.names = F, quote = F)

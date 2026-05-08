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

meta <- meta %>%
  mutate(Sample_Date = mdy(Sample_Collection_Date)) %>% # convert sampling dates into Date objects
  distinct(PopEnv, Sample_Date) %>%     # keep unique sampling dates 
  arrange(PopEnv, Sample_Date) %>% # Sort chronologically
  group_by(PopEnv) %>% # group by location so that we compute estimations independently within each location
  
  # previous sampling date
  mutate(Previous_Date = lag(Sample_Date)) %>% # create previous sampling date by shifting value down by one row
  
  ungroup()

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
    names_to = "Date",
    values_to = "ppt"
  ) %>% # convert to long format: each row = one day
  mutate(Date = ymd(Date)) # convert climate dates into Date objects

tmin_long <- tmin %>%
  pivot_longer(
    cols = matches("^2024"),
    names_to = "Date",
    values_to = "tmin"
  ) %>% # convert to long format: each row = one day
  mutate(Date = ymd(Date)) # convert climate dates into Date objects

tmean_long <- tmean %>%
  pivot_longer(
    cols = matches("^2024"),
    names_to = "Date",
    values_to = "tmean"
  ) %>% # convert to long format: each row = one day
  mutate(Date = ymd(Date)) # convert climate dates into Date objects

tmax_long <- tmax %>%
  pivot_longer(
    cols = matches("^2024"),
    names_to = "Date",
    values_to = "tmax"
  ) %>% # convert to long format: each row = one day
  mutate(Date = ymd(Date)) # convert climate dates into Date objects

sol_long <- sol %>%
  pivot_longer(
    cols = matches("^2024"),
    names_to = "Date",
    values_to = "sol"
  ) %>% # convert to long format: each row = one day
  mutate(Date = ymd(Date)) # convert climate dates into Date objects

vpdmax_long <- vpdmax %>%
  pivot_longer(
    cols = matches("^2024"),
    names_to = "Date",
    values_to = "vpdmax"
  ) %>% # convert to long format: each row = one day
  mutate(Date = ymd(Date)) # convert climate dates into Date objects

# -----------------------------
# 3. Calculate precipitation
#    between sampling dates
# -----------------------------
# the current sampling date IS included
# the previous sampling date is NOT included
# interval: (Previous_Date,Sample_Date]

climate_summary_ppt <- meta %>%
  left_join(ppt_long, by = "PopEnv", relationship = "many-to-many") %>% # join climate data to sampling intervals: match climate records with sampling intervals for the same location
  # relationship = "many-to-many" means that multiple rox of x/y match multiple rows of y/x
  # expected because we have several PopEnv (same sampling location) on different intervals (dates)
  filter(
    Date <= Sample_Date,
    Date > Previous_Date
  ) %>%
  group_by(PopEnv, Sample_Date, Previous_Date) %>%
  summarise(
    N_days       = n(),
    tot_ppt = sum(ppt, na.rm = TRUE),
    mean_ppt  = mean(ppt, na.rm = TRUE),
    .groups = "drop"
  )

climate_summary_tmin <- meta %>%
  left_join(tmin_long, by = "PopEnv", relationship = "many-to-many") %>% # join climate data to sampling intervals: match climate records with sampling intervals for the same location
  # relationship = "many-to-many" means that multiple rox of x/y match multiple rows of y/x
  # expected because we have several PopEnv (same sampling location) on different intervals (dates)
  filter(
    Date <= Sample_Date,
    Date > Previous_Date
  ) %>%
  group_by(PopEnv, Sample_Date, Previous_Date) %>%
  summarise(
    tmin = min(tmin, na.rm = TRUE),
    .groups = "drop"
  )

climate_summary_tmean <- meta %>%
  left_join(tmean_long, by = "PopEnv", relationship = "many-to-many") %>% # join climate data to sampling intervals: match climate records with sampling intervals for the same location
  # relationship = "many-to-many" means that multiple rox of x/y match multiple rows of y/x
  # expected because we have several PopEnv (same sampling location) on different intervals (dates)
  filter(
    Date <= Sample_Date,
    Date > Previous_Date
  ) %>%
  group_by(PopEnv, Sample_Date, Previous_Date) %>%
  summarise(
    tmean = mean(tmean, na.rm = TRUE),
    .groups = "drop"
  )

climate_summary_tmax <- meta %>%
  left_join(tmax_long, by = "PopEnv", relationship = "many-to-many") %>% # join climate data to sampling intervals: match climate records with sampling intervals for the same location
  # relationship = "many-to-many" means that multiple rox of x/y match multiple rows of y/x
  # expected because we have several PopEnv (same sampling location) on different intervals (dates)
  filter(
    Date <= Sample_Date,
    Date > Previous_Date
  ) %>%
  group_by(PopEnv, Sample_Date, Previous_Date) %>%
  summarise(
    tmax = max(tmax, na.rm = TRUE),
    .groups = "drop"
  )

climate_summary_sol <- meta %>%
  left_join(sol_long, by = "PopEnv", relationship = "many-to-many") %>% # join climate data to sampling intervals: match climate records with sampling intervals for the same location
  # relationship = "many-to-many" means that multiple rox of x/y match multiple rows of y/x
  # expected because we have several PopEnv (same sampling location) on different intervals (dates)
  filter(
    Date <= Sample_Date,
    Date > Previous_Date
  ) %>%
  group_by(PopEnv, Sample_Date, Previous_Date) %>%
  summarise(
    mean_sol = mean(sol, na.rm = TRUE),
    tot_sol = sum(sol, na.rm = TRUE),
    .groups = "drop"
  )

climate_summary_vpdmax <- meta %>%
  left_join(vpdmax_long, by = "PopEnv", relationship = "many-to-many") %>% # join climate data to sampling intervals: match climate records with sampling intervals for the same location
  # relationship = "many-to-many" means that multiple rox of x/y match multiple rows of y/x
  # expected because we have several PopEnv (same sampling location) on different intervals (dates)
  filter(
    Date <= Sample_Date,
    Date > Previous_Date
  ) %>%
  group_by(PopEnv, Sample_Date, Previous_Date) %>%
  summarise(
    mean_vpdmax = mean(vpdmax, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# 4. Join
# -----------------------------
climate_summary <- climate_summary_ppt %>%
  left_join(climate_summary_tmin %>% select(c("PopEnv", "Sample_Date", "Previous_Date"), tmin), by = c("PopEnv", "Sample_Date", "Previous_Date")) %>%
  left_join(climate_summary_tmean %>% select(c("PopEnv", "Sample_Date", "Previous_Date"), tmean), by = c("PopEnv", "Sample_Date", "Previous_Date")) %>%
  left_join(climate_summary_tmax %>% select(c("PopEnv", "Sample_Date", "Previous_Date"), tmax), by = c("PopEnv", "Sample_Date", "Previous_Date")) %>%
  left_join(climate_summary_sol %>% select(c("PopEnv", "Sample_Date", "Previous_Date"), c(mean_sol, tot_sol)), by = c("PopEnv", "Sample_Date", "Previous_Date")) %>%
  left_join(climate_summary_vpdmax %>% select(c("PopEnv", "Sample_Date", "Previous_Date"), mean_vpdmax), by = c("PopEnv", "Sample_Date", "Previous_Date")) 
  

# -----------------------------
# 5. Check
# -----------------------------
plot(climate_summary$Sample_Date, climate_summary$tmean, pch = 16)
plot(climate_summary$Sample_Date, climate_summary$mean_vpdmax, pch = 16)
plot(climate_summary$Sample_Date, climate_summary$tot_sol, pch = 16)


# -----------------------------
# 6. Export
# -----------------------------
write.table(climate_summary, 
            "/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/summaries/climate_intervals.csv", 
            sep = ",", col.names = T, row.names = F, quote = F)

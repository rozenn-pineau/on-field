#Goal: to extract and summarize the climate data 30 days before sampling times.

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

meta <- meta %>% mutate(Sample_Date = mdy(Sample_Collection_Date))

# -----------------------------
# 2. Get first sampling date
#    for each PopEnv
# -----------------------------

first_sampling <- meta %>%
  
  group_by(PopEnv) %>%
  summarise(First_Sample_Date = min(Sample_Date), .groups = "drop") %>% #find the minimum date within the group
  mutate(Start_Date = First_Sample_Date - days(30)) #remove 30 days (days() create 30 days in date format)

# -----------------------------
# 3. Load climate datasets
# -----------------------------
ppt <- fread("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/extracted/ppt_April_October_2024.csv")
tmin <- fread("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/extracted/tmin_April_October_2024.csv")
tmean <- fread("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/extracted/tmean_April_October_2024.csv")
tmax <- fread("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/extracted/tmax_April_October_2024.csv")
sol <- fread("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/extracted/sol_April_October_2024.csv")
vpdmax <- fread("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/extracted/vpdmax_April_October_2024.csv")

# -----------------------------
# 4. convert from wide to long format
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
# 5. Extract climate during
#    the 30 days before sampling
# -----------------------------

ppt_30d <- first_sampling %>%
  left_join(ppt_long, by = "PopEnv") %>%
  filter(Date >= Start_Date, Date < First_Sample_Date) %>% #interval: [First Sample−30, First Sample)
  group_by(PopEnv, First_Sample_Date) %>%
  summarise(tot_ppt_30d_prior = sum(ppt, na.rm = TRUE),
    mean_ppt_30d_prior  = mean(ppt, na.rm = TRUE), .groups = "drop")

tmin_30d <- first_sampling %>%
  left_join(tmin_long, by = "PopEnv") %>%
  filter(Date >= Start_Date, Date < First_Sample_Date) %>% #interval: [First Sample−30, First Sample)
  group_by(PopEnv, First_Sample_Date) %>%
  summarise(tmin_30d_prior = min(tmin, na.rm = TRUE),.groups = "drop")

tmean_30d <- first_sampling %>%
  left_join(tmean_long, by = "PopEnv") %>%
  filter(Date >= Start_Date, Date < First_Sample_Date) %>% #interval: [First Sample−30, First Sample)
  group_by(PopEnv, First_Sample_Date) %>%
  summarise(tmean_30d_prior = mean(tmean, na.rm = TRUE),.groups = "drop")

tmax_30d <- first_sampling %>%
  left_join(tmax_long, by = "PopEnv") %>%
  filter(Date >= Start_Date, Date < First_Sample_Date) %>% #interval: [First Sample−30, First Sample)
  group_by(PopEnv, First_Sample_Date) %>%
  summarise(tmax_30d_prior = max(tmax, na.rm = TRUE),.groups = "drop")

sol_30d <- first_sampling %>%
  left_join(sol_long, by = "PopEnv") %>%
  filter(Date >= Start_Date, Date < First_Sample_Date) %>% #interval: [First Sample−30, First Sample)
  group_by(PopEnv, First_Sample_Date) %>%
  summarise(tot_sol_30d_prior = sum(sol, na.rm = TRUE),
            mean_sol_30d_prior = mean(sol, na.rm = TRUE), .groups = "drop")

vpdmax_30d <- first_sampling %>%
  left_join(vpdmax_long, by = "PopEnv") %>%
  filter(Date >= Start_Date, Date < First_Sample_Date) %>% #interval: [First Sample−30, First Sample)
  group_by(PopEnv, First_Sample_Date) %>%
  summarise(vpdmax_30d_prior = mean(vpdmax, na.rm = TRUE),.groups = "drop")

# -----------------------------
# 6. Join
# -----------------------------
climate_summary_30d_prior <- ppt_30d %>%
  left_join(tmin_30d %>% select(c("PopEnv", "First_Sample_Date"), tmin_30d_prior), by = c("PopEnv", "First_Sample_Date")) %>%
  left_join(tmean_30d %>% select(c("PopEnv", "First_Sample_Date"), tmean_30d_prior), by = c("PopEnv", "First_Sample_Date")) %>%
  left_join(tmax_30d %>% select(c("PopEnv", "First_Sample_Date"), tmax_30d_prior), by = c("PopEnv", "First_Sample_Date")) %>%
  left_join(sol_30d %>% select(c("PopEnv", "First_Sample_Date"), c(tot_sol_30d_prior, mean_sol_30d_prior)), by = c("PopEnv", "First_Sample_Date")) %>%
  left_join(vpdmax_30d %>% select(c("PopEnv", "First_Sample_Date"), vpdmax_30d_prior), by = c("PopEnv", "First_Sample_Date")) 

# -----------------------------
# 7. Check
# -----------------------------
plot(climate_summary_30d_prior$First_Sample_Date, climate_summary_30d_prior$tmin_30d_prior, pch = 16)
plot(climate_summary_30d_prior$First_Sample_Date, climate_summary_30d_prior$tot_ppt_30d_prior, pch = 16)
plot(climate_summary_30d_prior$First_Sample_Date, climate_summary_30d_prior$tot_sol_30d_prior, pch = 16)

# -----------------------------
# 6. Export
# -----------------------------
write.table(climate_summary_30d_prior, 
            "/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/summaries/climate_30d_prior_to_first_sampling.csv", 
            sep = ",", col.names = T, row.names = F, quote = F)


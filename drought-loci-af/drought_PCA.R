#Goal: to explore what the data looks like, and to plot a few first analyses through space and time

#Clear environment
rm(list = ls())

#Libraries
library(tidyverse)
library(ggplot2)
library(vegan) #for pRDA
library(data.table)
library(car)
library(tidyverse)
library(broom)
library(lme4) #for linear mixed models
library(lmerTest) #for linear mixed models
library(broom.mixed) #for linear mixed models
library(dplyr)
library(tidyr)
library(purrr)


# -----------------------------
# Define constants
# -----------------------------
n_pools <- 87
#Outliers to remove
outliers <- c("JK-JM-Ohi-pl33_S55", "JK-JM-Ohi-pl36_S58")
#I ran the PCA a first time and those two pools were very much different from the others, we remove them from then on

#exclude Ontario samples for now (no climate variables)
ontario <- c("Ontario", "Ontario 1")

# -----------------------------
# Load pool allele frequencies
# -----------------------------
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/drought_loci")
loci <- fread("drought_loci_AF_matrix_2024_89pools.csv")
loci_clean <- loci[, !..outliers] #for data.table, which read expressions rather than column numbers, 
#".." this will read all columns except the ones named in drop_cols

# -----------------------------
# Load metadata
# -----------------------------
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/metadata")
metadata <- fread("Sample_info.csv")
metadata <- as.data.frame(metadata)
head(metadata)
metadata_clean <- metadata[metadata$Sample!=outliers,]
metadata_clean <- metadata_clean[!metadata_clean$Population %in% ontario,] 

# -----------------------------
# Load climate data
# -----------------------------
clim_single <- fread("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/summaries/climate_time_points_w_PCA.csv")
head(clim_single)
clim_single <- as.data.frame(clim_single)
clim_single_clean <- clim_single[clim_single$Sample!=outliers, ]
clim_30d <- fread("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/summaries/climate_30d_prior_to_each_sampling_w_PCA.csv")
head(clim_30d)
clim_30d <- as.data.frame(clim_30d)
clim_30d_clean <- clim_30d[!clim_30d$Sample %in% outliers, ]

# -----------------------------
# Long to wide format
# -----------------------------
#long to wide format for loci to be able to join with metadata
loci_long <- loci_clean %>%
  pivot_longer(
    cols = -c(snp_id, chrom, pos, pos2, ref),   # everything EXCEPT the SNP info columns
    names_to = "pool_name",             
    values_to = "pool_AF"              
  ) # 893 SNPs x 87 pools is 77691 rows

# -----------------------------
# Merge metadata
# -----------------------------
merged <- loci_long %>% left_join(metadata_clean, by = c("pool_name" = "Sample")) # 77691 - good
head(merged)

# -----------------------------
# Merge climate data
# -----------------------------
  #metadata
metadata_w_climate_single_pt  <- metadata_clean %>% left_join(clim_single_clean, by = "Sample") #for PCA coloring
  #AF
merged_w_climate_single_pt <- merged %>% left_join(clim_single_clean, by = c("pool_name" = "Sample")) # 77691 rows - good
merged_w_climate30d <- merged %>% left_join(clim_30d_clean, by = c("pool_name" = "Sample")) # 77691 rows - good


# --------------------------------------
# Prepare the data for PCA
# --------------------------------------

#prepare the matrix
freq_mat <- loci_clean[,c(6:92)] %>% t()   # 87 pools x 893 SNPs

#check missingness
na_per_snp  <- colSums(is.na(freq_mat))
na_per_pool <- rowSums(is.na(freq_mat))
summary(na_per_snp)

#Drop SNPs with more than 20% missing data
freq_mat_clean <- freq_mat[,-which(na_per_snp/n_pools >= 0.2)] #89 pools x 877 SNPs

# Impute: replace each NA with its SNP's mean across pools
freq_mat_imputed <- apply(freq_mat_clean, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})


# --------------------------------------
# PCA
# --------------------------------------
pca <- prcomp(freq_mat_imputed, center = TRUE, scale. = FALSE)

#How many PCs matter?
screeplot(pca, npcs = 10) #mostly PC1 and PC2

# Percent variance explained by each PC
var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)

# Match metadata to PCA row order
meta_ord <- metadata_w_climate_single_pt[match(rownames(pca$x), metadata_w_climate_single_pt$Sample), ]

# --------------------------------------
# PCA plot
# --------------------------------------
df <- data.frame(pca$x[, 1:4], meta_ord) #make a new df with the 4 first PCs and the metadata

ggplot(df, aes(PC1, PC2, color = Environment)) +   # try location, lat, precipitation, etc.
  geom_point(size = 2) +
  coord_equal() +
  labs(x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2]))


#test for longitudinal effect on PC1 - simple regression
summary(lm(PC1 ~ long, data = df))
summary(lm(PC1 ~ Population, data = df))
summary(lm(PC1 ~ ppt, data = df))
summary(lm(PC1 ~ tmean, data = df))
summary(lm(PC1 ~ tmax, data = df))
summary(lm(PC1 ~ tmin, data = df))
summary(lm(PC1 ~ sol, data = df))
summary(lm(PC1 ~ vpdmax, data = df))
summary(lm(PC1 ~ PC1_clim, data = df))
summary(lm(PC1 ~ PC2_clim, data = df))
summary(lm(PC1 ~ Environment, data = df))
summary(lm(PC1 ~ CornSoy, data = df))

# --------------------------------------
# Change through time
# --------------------------------------
ggplot(merged, aes(x = Month, y = pool_AF)) +
  geom_line(aes(group = snp_id), alpha = 0.05, color = "grey40") +
  stat_summary(fun = mean, geom = "line", color = "firebrick", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", color = "firebrick", size = 2) +
  facet_wrap(~ Population) +
  labs(x = "Sampling month", y = "Allele frequency") +
  theme_bw()


merged_w_climate_single_pt %>%
  group_by(Environment, Population, Month) %>%
  summarise(mean_freq = mean(pool_AF, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(Month, mean_freq, color = Population, group = Population)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Environment, ncol = 2) +
  labs(x = "Sampling month", y = "Mean allele frequency") +
  theme_bw()

# --------------------------------------
# per population mean AF versus climate change
# --------------------------------------

#calculate the correlation between change in climate and change in mean allele frequency per population
#for each time point where we have a mean AF per population, what is the climate for this single point?

#graph example - AF as a function of time
merged_w_climate_single_pt %>%
  filter(Population == "DPAC") %>%
  group_by(Environment, Population, Month) %>%
  summarise(mean_freq = mean(pool_AF, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(Month, mean_freq, color = Population, group = Population)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Environment, ncol = 2) +
  labs(x = "Sampling month", y = "Mean allele frequency") +
  theme_bw()


#climate as a function of time
merged_w_climate_single_pt %>%
  filter(Population == "DPAC") %>%
  group_by(Environment, Population, Month) %>%
  ggplot(aes(Month, tmean, color = Population, group = Population)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Environment, ncol = 2) +
  labs(x = "Sampling month", y = "Mean temperature") +
  theme_bw()


#AF as a function of climate
merged_w_climate_single_pt %>%
  filter(Population == "DPAC") %>%
  group_by(Environment, Population, Month, tmean) %>%
  summarise(mean_freq = mean(pool_AF, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(tmean, mean_freq, color = Population, group = Population)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Environment, ncol = 2) +
  labs(x = "Mean temperature", y = "Mean allele frequency") +
  theme_bw()

# --------------------------------------
# single climate point
# --------------------------------------

#Models of AF versus climate by population, single time point climate for now (use 30 days prior to each sampling next)

clim_vars <- c("ppt", "tmin", "tmean", "tmax", "sol", "vpdmax", "PC1_clim")

#Calculate the mean AF per sample/time point (keep climate values alongside)
af_summary <- merged_w_climate_single_pt %>%
  group_by(Population, Environment, Month, Year) %>%
  summarise(mean_freq = mean(pool_AF, na.rm = TRUE),
    across(all_of(clim_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

#Long format: one row per climate variable
af_long <- af_summary %>% pivot_longer(all_of(clim_vars), names_to = "climate_var", values_to = "clim_value")

#Correlation test within each Population x Environment x climate variable
#Exclude Ontario (no climate variable)

cor_results <- af_long %>%
  group_by(Population, Environment, climate_var) %>%
  filter(n() >= 3) %>%                      #only test of more than 3 time points
  summarise(n = n(),
            test = list(cor.test(clim_value, mean_freq, method = "spearman")), .groups = "drop") %>%
  mutate(tidy = map(test, tidy)) %>%
  unnest(tidy) %>% 
  select(Population, Environment, climate_var, n, r = estimate,p.value) %>% #, conf.low, conf.high
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%  # multiple testing
  arrange(p_adj)

cor_results$p_adj < 0.05
#correlation analysis on single time point for climate and mean allele frequency did not give out anything significant 

#Try the same analysis on the climate data 30d prior
# --------------------------------------
# climate 30d prior
# --------------------------------------
clim_vars <- c("tot_ppt_30d_prior", "mean_ppt_30d_prior", "tmin_30d_prior", "tmean_30d_prior", 
               "tmax_30d_prior", "tot_sol_30d_prior", "mean_sol_30d_prior",  "vpdmax_30d_prior", "PC1_clim" )

#Calculate the mean AF per sample/time point (keep climate values alongside)
af_summary <- merged_w_climate30d %>%
  group_by(Population, Environment, Month, Year) %>%
  summarise(mean_freq = mean(pool_AF, na.rm = TRUE),
            across(all_of(clim_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

#Long format: one row per climate variable
af_long <- af_summary %>% pivot_longer(all_of(clim_vars), names_to = "climate_var", values_to = "clim_value")

#Correlation test within each Population x Environment x climate variable
cor_results <- af_long %>%
  group_by(Population, Environment, climate_var) %>%
  filter(n() >= 3) %>%                      #only test of more than 3 time points
  summarise(n = n(),
            test = list(cor.test(clim_value, mean_freq, method = "spearman")), .groups = "drop") %>%
  mutate(tidy = map(test, tidy)) %>%
  unnest(tidy) %>% 
  select(Population, Environment, climate_var, n, r = estimate,p.value) %>% #, conf.low, conf.high
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%  # multiple testing
  arrange(p_adj)

cor_results$p_adj < 0.05

#This version of climate data also did not give any significant results


# --------------------------------------
# Linear mixed model with random and fixed effect on all loci at once
# --------------------------------------

#Explanation for this code loop
#Groups the data by Environment and climate variable. It then run a mixed effect model
#on the mean frequency, testing for the effect of climate and controlling for the random effect
#of population.
#It extracts the statistics and p-value for the model. 
#To note, the order in whixh libraries are called for lme4 and lmerTest matters, 
#lme4 alone will not give p-values (ambiguous degrees of freedom). 
# lmerTest package fixes this by wrapping lmer() and adding approximate degrees of freedom (Satterthwaite).


# Trying when ag and nat are tested separately
mixed_results <- af_long %>%
  drop_na(clim_value, mean_freq) %>%
  group_by(Environment, climate_var) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lmer(mean_freq ~ scale(clim_value) + (1 | Population), data = .x)),
    tidy  = map(model, ~ tidy(.x, effects = "fixed"))
  ) %>%
  unnest(tidy) %>%
  filter(term == "scale(clim_value)") %>%
  select(Environment, climate_var, estimate, std.error, statistic, any_of("p.value")) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))

mixed_results

# Trying when ag and nat are tested together
mixed_results <- af_long %>%
  drop_na(clim_value, mean_freq) %>%
  group_by(climate_var) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lmer(mean_freq ~ scale(clim_value) + (1 | Population), data = .x)),
    tidy  = map(model, ~ tidy(.x, effects = "fixed"))
  ) %>%
  unnest(tidy) %>%
  filter(term == "scale(clim_value)") %>%
  select(climate_var, estimate, std.error, statistic, any_of("p.value")) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))

mixed_results

#I did not detect any significant relationship but it might be because we do not have enough power

# --------------------------------------
# SNP by SNP analysis
# --------------------------------------
#why: maybe if we do not detect anything, is because there are loci that are not responding 
#to env changes (878 loci may blurr the signal?)
#--> downsample to significantly changing loci?

#dataset with 878 loci and climate data
head(merged_w_climate_single_pt)

#identify loci between populations that change significantly

#example with one population and one SNP
test <- merged_w_climate_single_pt %>% filter(snp_id == "1", Population == "DPAC", Environment == "ag")
summary(lm(test$Month ~ test$pool_AF))


fit_one <- function(d) {
  # skip series that are too short or have no variation
  if (nrow(d) < 4 || n_distinct(d$Month) < 3 || var(d$pool_AF) == 0) return(NULL)
  
  m <- lm(pool_AF ~ Month, data = d)
  tidy(m) %>%
    filter(term == "Month") %>%
    transmute(slope = estimate, se = std.error, t = statistic,
              p = p.value, n_obs = nrow(d))
}

results <- merged_w_climate_single_pt %>%
  group_by(snp_id, Population, Environment) %>%
  nest() %>%
  mutate(res = map(data, possibly(fit_one, otherwise = NULL))) %>%
  select(-data) %>%
  unnest(res) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p, method = "BH"))

head(results)
which(results$p_adj < 0.05)
table(results[which(results$p < 0.05),])


# counts per SNP
snp_counts <- table(results$snp_id[which(results$p < 0.05)])

# how many SNPs appear 1 time, 2 times, 3 times, ...
table(snp_counts)

# the SNPs that appear most often
sort(snp_counts, decreasing = TRUE) %>% head(10)

#SNPs 205 236 268 392 716 significantly change in three different populations (before pval correction)

#plot and test against climate in those locations?
#is this different than what you would expect by chance? (do the same analyses on a set of random samples)

#To try 
#Per-SNP models with FDR correction: fit AF ~ tmean + Population for each SNP (or (1 | Population)), 
# extract the tmean p-values, and apply p.adjust(method = "fdr"). This tells you which SNPs respond.

#Partial RDA , rda(freq_imp ~ tmean + Condition(Population)) 
#with restricted permutations, tests all SNPs jointly without the pseudoreplication problem.







#identify loci linked with climate using two different approaches
#single models with drought SNPs - consistently changing candidtae loci within the pool of candidate loci
#GEA with climate data

#Does it make sense to look at the drought PGS and potential association with climate?
#multiple regression
#--> what model structure to try?
#try one model per SNP as in Battlay 2026

#PGS score on pooled Seq data makes sense











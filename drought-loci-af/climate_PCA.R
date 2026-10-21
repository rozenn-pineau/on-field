#Goal: to explore what the data looks like, and to plot a few first analyses through space and time

#Libraries
library(tidyverse)
library(ggplot2)
library(vegan) #for pRDA
library(data.table)
library(car)
library(permute)
library(Hmisc) #for p-values of correlation matrix
library(corrplot) #for correlation plots


#Clear environment
rm(list = ls())


# -----------------------------
# Load pool allele frequencies
# -----------------------------
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/drought_loci")
loci <- fread("drought_loci_AF_matrix_2024_89pools.csv", sep = ",")
summary(loci)

# -----------------------------
# Load climate data
# -----------------------------
clim_single <- fread("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/summaries/climate_time_points.csv")
summary(clim_single)
clim_30d <- fread("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/summaries/climate_30d_prior_to_each_sampling.csv")
summary(clim_30d)

# -----------------------------
# Define constants
# -----------------------------
n_pools <- 89

# --------------------------------------
# What does the climate data look like? --> PCA on climate data
# --------------------------------------

# --------------------------------------
# Single time points dataset
# --------------------------------------

#prepare the matrix
climate_mat <- clim_single[,19:24]   # 96 pools x 893 SNPs - need to filter out some pools
rownames(climate_mat) <- clim_single$Sample
climate_mat_89pools <- climate_mat[which(rownames(climate_mat) %in% colnames(loci)[6:94]),] # 89 pools x 6 climate variables
rownames(climate_mat_89pools) <- rownames(climate_mat)[which(rownames(climate_mat) %in% colnames(loci)[6:94])]

#check that I kept the correct names
rownames(climate_mat_89pools) %in% colnames(loci) #good

#remove empty rows - Ontario and 2025 pools
climate_mat_89pools_clean <- climate_mat_89pools[rowSums(is.na(climate_mat_89pools))<=1,] #81 pools x 6 climate variables
rownames(climate_mat_89pools_clean) <- rownames(climate_mat_89pools)[rowSums(is.na(climate_mat_89pools))<=1]

# Impute: replace each NA with the climate mean across pools
climate_mat_89pools_clean_imputed <- apply(climate_mat_89pools_clean, 2, function(x) { #summarize mean across colums
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})

#add pool names
rownames(climate_mat_89pools_clean_imputed) <- rownames(climate_mat_89pools_clean)

#clean environment
rm(climate_mat_89pools, climate_mat, climate_mat_89pools_clean)

# --------------------------------------
# PCA
# --------------------------------------
pca <- prcomp(climate_mat_89pools_clean_imputed, center = TRUE, scale. = TRUE)

#How many PCs matter?
screeplot(pca, npcs = 10) #mostly PC1 and PC2

# Percent variance explained by each PC
var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)

# --------------------------------------
# PCA plot
# --------------------------------------
df <- data.frame(pca$x[, 1:4], climate_mat_89pools_clean_imputed) #add initial values for coloring

ggplot(df, aes(PC1, PC2, color = ppt)) +   
  geom_point(size = 2) +
  coord_equal() +
  labs(x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2]), 
       title = "*Day of sampling* climate PCA")

#export PC1 and PC2
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/summaries")
#write.table(df, "climate_time_points_w_PCA.csv", sep = ",", col.names = T, row.names = T, quote = F)
#Make the PCA for all three climate datasets, export PC1 and PC2 to test on associations with genotypes

# --------------------------------------
# Correlation between variables
# --------------------------------------
#prepare the matrix
clim <- df %>% select(ppt, tmin, tmean, tmax, sol, vpdmax)

#Correlation matrix (Spearman is robust to skewed variables like ppt)
cor_mat <- cor(clim, method = "spearman")
round(cor_mat, 2)

#Matrix with p-values
res <- rcorr(as.matrix(clim), type = "spearman")
res$r   # correlation coefficients
res$P   # p-values

#Visualise
corrplot(cor_mat, method = "color", type = "upper",
         addCoef.col = "black", tl.col = "black",
         p.mat = res$P, sig.level = 0.05, insig = "blank")

rm(df, cor_mat, pca, res)

# --------------------------------------
# *30 days prior to each sampling* dataset
# --------------------------------------
#prepare the matrix
climate_mat <- clim_30d[,c(1,7:14)]   # 90 pools x 8 climate variables 
#pools filtered out in this dataset: Ontario (4) and Tazwell corn mispelled (2)

#filter out pools that are filtered out in loci (89 poos)
# For every element in x, is it present anywhere inside y
climate_mat_84pools <- climate_mat[which(climate_mat$Sample %in% colnames(loci)[6:94]),] # 84 pools x 8 climate variables

#remove empty rows - Ontario and 2025 pools
climate_mat_84pools_clean <- climate_mat_84pools[rowSums(is.na(climate_mat_84pools))<=1,] # 81 pools x 8 climate variables

# Impute: replace each NA with the climate mean across pools
climate_mat_84pools_clean_imputed <- apply(climate_mat_84pools_clean[,c(2:9)], 2, function(x) { #summarize mean across colums
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})

#add pool names
climate_mat_84pools_clean_imputed_rownames <- data.frame(Sample = climate_mat_84pools_clean$Sample, 
                                                         climate_mat_84pools_clean_imputed) 

#clean environment
rm(climate_mat_84pools, climate_mat, climate_mat_84pools_clean, climate_mat_84pools_clean)

# --------------------------------------
# PCA
# --------------------------------------
pca <- prcomp(climate_mat_84pools_clean_imputed_rownames[,c(2:9)], center = TRUE, scale. = TRUE)

#How many PCs matter?
screeplot(pca, npcs = 10) #mostly PC1 and PC2

# Percent variance explained by each PC
var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)

# --------------------------------------
# PCA plot
# --------------------------------------
df <- data.frame(Sample = climate_mat_84pools_clean_imputed_rownames$Sample, 
                 pca$x[, 1:4], 
                 climate_mat_84pools_clean_imputed_rownames[,2:9]) #add initial values for coloring

ggplot(df, aes(PC1, PC2, color = tmean_30d_prior)) +   
  geom_point(size = 2) +
  coord_equal() +
  labs(x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2]), 
       title = "*30 days before sampling* climate PCA")

#export PC1 and PC2
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/PRISM/summaries")
#write.table(df, "climate_30d_prior_to_each_sampling_w_PCA.csv", sep = ",", col.names = T, row.names = F, quote = F)
#Make the PCA for all three climate datasets, export PC1 and PC2 to test on associations with genotypes

# --------------------------------------
# Correlation between variables
# --------------------------------------
#prepare the matrix
clim <- df %>% select(tot_ppt_30d_prior, mean_ppt_30d_prior, tmin_30d_prior, tmean_30d_prior, tmax_30d_prior,   
                      tot_sol_30d_prior, mean_sol_30d_prior, vpdmax_30d_prior)

#Correlation matrix (Spearman is robust to skewed variables like ppt)
cor_mat <- cor(clim, method = "spearman")
round(cor_mat, 2)

#Matrix with p-values
res <- rcorr(as.matrix(clim), type = "spearman")
res$r   # correlation coefficients
res$P   # p-values

#Visualise
corrplot(cor_mat, method = "color", type = "upper",
         addCoef.col = "black", tl.col = "black",
         p.mat = res$P, sig.level = 0.05, insig = "blank")


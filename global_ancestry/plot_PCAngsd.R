#!/usr/bin/env Rscript
# working directory is hilo/
# plots PCA results
library(dplyr)
library(tidyr)
library(ggplot2)

# load variables from snakemake
# get covariance matrix estimated by PCAngsd
cov_file = snakemake@input[["cov"]]
# cov_file = "global_ancestry/results/PCA/HILO_MAIZE55/whole_genome.cov"

# get output plot filenames
png_pca = snakemake@output[["png_pca"]]
# png_pca = "global_ancestry/plots/pca.png"
png_pc34 = snakemake@output[["png_pc34"]]
# png_pc34 = "global_ancestry/plots/pc34.png"
png_pc56 = snakemake@output[["png_pc56"]]
# png_pc56 = "global_ancestry/plots/pc56.png"

# get colors for plot
# source("colors.R")
source(snakemake@input[["colors"]])

# get sample metadata, including estimated sequence coverage
load(snakemake@input[["meta"]]) 
# load("samples/HILO_MAIZE55_meta.RData")
# note: metadata is in the same order as samples input to PCAngsd


# get PCA data
cov_data <- read.table(cov_file, header = F, stringsAsFactors = F)
# PC's are each column of dataframe pca:
# i.e. PC1 is V1, PC2 is V2 etc.
pca <- eigen(cov_data) # take PCA of covariance matrix
m = 10 # make small dataframe w/ only first m eigenvectors
# eigenvectors
pca_small <- data.frame(pca$vectors[ , 1:m])
colnames(pca_small) = paste0("PC", 1:m)
# rounded eigen values
PC_var_explained = 100*pca$values/sum(pca$values)

# join ids and first PCs of covariance PCA data by position (CAUTION - bam list order and admix results MUST MATCH!)
d <- bind_cols(meta, pca_small)  %>%
  arrange(., popN) %>%
  arrange(., zea) %>%
  arrange(., symp_allo) %>%
  arrange(., group, ELEVATION)

# plot first 2 PC's
p12 = d %>%
  ggplot(., aes(PC1, PC2)) + 
  xlab(paste0("PC1 (", round(PC_var_explained[1], 2), "%)")) +
  ylab(paste0("PC2 (", round(PC_var_explained[2], 2), "%)")) +
  theme_classic() +
  geom_point(aes(color = group, 
                 shape = group,
                 size = log10(est_coverage)), 
             alpha = 0.5) +
  scale_color_manual(values = col_group_zea, labels = zea_group_labels) +
  scale_shape_manual(values = shape_group_zea, labels = zea_group_labels) +
  scale_size_continuous(breaks = -2:3, labels = paste0(10^(-2:3), "x")) +
  labs(color = "Zea", shape = "Zea", size = "Coverage")
#p12
ggsave(png_pca, 
       plot = p12, 
       device = "png", 
       width = 5.4, height = 3.5, units = "in",
       dpi = 200)

p34 = d %>%
  ggplot(., aes(PC3, PC4)) + 
  xlab(paste0("PC3 (", round(PC_var_explained[3], 2), "%)")) +
  ylab(paste0("PC4 (", round(PC_var_explained[4], 2), "%)")) +
  theme_classic() +
  geom_point(aes(color = group, 
                 shape = group,
                 size = log10(est_coverage)), 
             alpha = 0.5) +
  scale_color_manual(values = col_group_zea, labels = zea_group_labels) +
  scale_shape_manual(values = shape_group_zea, labels = zea_group_labels) +
  scale_size_continuous(breaks = -2:3, labels = paste0(10^(-2:3), "x")) +
  labs(color = "Zea", shape = "Zea", size = "Coverage")
#p34
ggsave(png_pc34, 
       plot = p34, 
       device = "png", 
       width = 5.4, height = 3.5, units = "in",
       dpi = 200)
p56 = d %>%
  ggplot(., aes(PC5, PC6)) + 
  xlab(paste0("PC5 (", round(PC_var_explained[5], 2), "%)")) +
  ylab(paste0("PC6 (", round(PC_var_explained[6], 2), "%)")) +
  theme_classic() +
  geom_point(aes(color = group, 
                 shape = group,
                 size = log10(est_coverage)), 
             alpha = 0.5) +
  scale_color_manual(values = col_group_zea, labels = zea_group_labels) +
  scale_shape_manual(values = shape_group_zea, labels = zea_group_labels) +
  scale_size_continuous(breaks = -2:3, labels = paste0(10^(-2:3), "x")) +
  labs(color = "Zea", shape = "Zea", size = "Coverage")
p56
ggsave(png_pc56, 
       plot = p56, 
       device = "png", 
       width = 5.4, height = 3.5, units = "in",
       dpi = 200)
#!/usr/bin/env Rscript
# working directory is hilo/
# plots PCA results for putative inversion at mhl1
library(dplyr)
library(tidyr)
library(ggplot2)

# load variables from snakemake
# get covariance matrix estimated by PCAngsd
cov_file = snakemake@input[["cov"]]
# cov_file = "mhl1_inv/results/PCA/HILO_MAIZE55/Ne10000_yesBoot/mhl1_inv.cov"

# get output plot filenames
png_pca = snakemake@output[["png_pca"]]
# png_pca = "mhl1_inv/plots/HILO_MAIZE55/Ne10000_yesBoot/mhl1_inv_pca.png"
png_pca_lzw = snakemake@output[["png_pca_lzw"]]
# png_pca_lzw = "../hilo_manuscript/figures_supp/Ne10000_yesBoot_mhl1_inv_pca.tif"

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

# update labels for plot legend:
zea_group_labels["allopatric_maize"] <- "Reference maize"
zea_group_labels["allopatric_mexicana"] <- "Reference mexicana"
zea_group_labels["parviglumis"] <- "Reference parviglumis"

# plot first 2 PC's
p12 = d %>%
  ggplot(., aes(PC1, PC2)) + 
  xlab(paste0("PC1 (", round(PC_var_explained[1], 2), "%)")) +
  ylab(paste0("PC2 (", round(PC_var_explained[2], 2), "%)")) +
  theme_classic() +
  geom_point(aes(color = group,
                 alpha = group,
                 shape = group
                 )) +
  scale_alpha_manual(values = alphas_group_zea, labels = zea_group_labels) +
  scale_color_manual(values = col_group_zea, labels = zea_group_labels) +
  scale_shape_manual(values = shape_group_zea, labels = zea_group_labels) +
  labs(color = "Zea", shape = "Zea", alpha = "Zea")
# p12
ggsave(png_pca, 
       plot = p12, 
       device = "png", 
       width = 6, height = 3.5, units = "in",
       dpi = 300)

ggsave(png_pca_lzw, 
       plot = p12, 
       device = "tiff", 
       width = 6, height = 3.5, units = "in",
       dpi = 300, 
       compression = "lzw", type = "cairo")
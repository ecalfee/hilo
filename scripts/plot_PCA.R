# plots PCA results
library(scales)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)# for color palette

# standard colors and labels for plots
blues = brewer.pal(n = 9, name = "Blues")[c(6,8)]
yellows = brewer.pal(n = 9, name = "YlOrBr")[c(4,6)]
colors_maize2mex = c(yellows, blues)

# get PCA data
pca_dir <- "../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedPCA/PCA"
pca <- read.table(file.path(pca_dir, "whole_genome.cov"),
                       header = F, stringsAsFactors = F)
# PC's are each column of dataframe pca:
# i.e. PC1 is V1, PC2 is V2 etc.

# get metadata for individuals included in NGSadmix analysis
pass1_allo4Low <- read.table("../data/pass1_allo4Low_ids.txt", stringsAsFactors = F, 
                             header = T, sep = "\t")
# join bams and first 10 PCs of covariance PCA data by position (CAUTION - bam list order and admix results MUST MATCH!)
d <- bind_cols(pass1_allo4Low, pca[ ,1:10])  %>%
  arrange(., popN) %>%
  arrange(., zea) %>%
  arrange(., symp_allo) %>%
  mutate(., group = paste(symp_allo, zea, sep = "_"))

# plot PCA


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
colors_alphabetical = colors_maize2mex[c(1,4,2,3)] # allo maize, allo mex, symp maize, symp mex

# get PCA data
cov_dir <- "../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedPCA/PCA"
cov_data <- read.table(file.path(cov_dir, "whole_genome.cov"),
                       header = F, stringsAsFactors = F)
# PC's are each column of dataframe pca:
# i.e. PC1 is V1, PC2 is V2 etc.
pca <- eigen(cov_data) # take PCA of covariance matrix
n = 4 # make small dataframe w/ only first n eigenvectors
pca_small <- data.frame(pca$vectors[ , 1:n])
colnames(pca_small) = paste0("PC", 1:n)
  
# get metadata for individuals included in NGSadmix analysis
pass1_allo4Low <- read.table("../data/pass1_allo4Low_ids.txt", stringsAsFactors = F, 
                             header = T, sep = "\t")
# join bams and first 10 PCs of covariance PCA data by position (CAUTION - bam list order and admix results MUST MATCH!)
d <- bind_cols(pass1_allo4Low, pca_small)  %>%
  #arrange(., popN) %>%
  #arrange(., zea) %>%
  #arrange(., symp_allo) %>%
  mutate(., group = paste(symp_allo, zea, sep = "_"))

# rounded eigen values
PC_var_explained = round(pca$values, 2)

# plot first PC's
# PC1 and 2, all no filter
p12_coverage = d %>%
  ggplot(., aes(PC1, PC2, alpha = .5)) + 
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC1 (", PC_var_explained[2], "%)")) +
  ggtitle("PC1 and PC2 (all ind's); SNPs spaced >.01cM") +
  geom_point(aes(color = group, size = est_coverage)) +
  scale_colour_manual(values = colors_alphabetical)
ggsave("../plots/PCA_12_by_coverage.png", plot = p12_coverage, device = png(), 
       width = 12, height = 8, units = "in",
       dpi = 200)
# PC1 and 2, excluding super low coverae ind's
p12 = d %>%
  filter(., est_coverage >= .05) %>%
  ggplot(., aes(PC1, PC2)) + 
  ggtitle("PC1 and PC2 (excl. low coverage ind's < .05x); SNPs spaced >.01cM") +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC1 (", PC_var_explained[2], "%)")) +
  geom_point(aes(color = group)) +
  scale_colour_manual(values = colors_alphabetical)
plot(p12)
ggsave("../plots/PCA_12.png", plot = p12, device = png(), 
       width = 12, height = 8, units = "in",
       dpi = 200)
# PC3 and PC4, excluding low coverage ind's
p34 = d %>%
  filter(., est_coverage >= .05) %>%
  ggplot(., aes(PC3, PC4)) + 
  ggtitle("PC3 and PC4 (excl. low coverage ind's < .05x); SNPs spaced >.01cM") +
  xlab(paste0("PC3 (", PC_var_explained[3], "%)")) +
  ylab(paste0("PC4 (", PC_var_explained[4], "%)")) +
  geom_point(aes(color = group)) +
  scale_colour_manual(values = colors_alphabetical)
plot(p34)
ggsave("../plots/PCA_34.png", plot = p34, device = png(), 
       width = 12, height = 8, units = "in",
       dpi = 200)

# Use vars() to supply faceting variables:
p12 + facet_wrap(~LOCALITY)
ggsave("../plots/PCA_facet_by_population.png", 
       plot = p12 + facet_wrap(~LOCALITY), device = png(), 
       width = 12, height = 8, units = "in",
       dpi = 200)


# repeat same analyses for SNP data pruned every 1000 SNPs
cov_data_1000 <- read.table("../data/geno_lik/merged_pass1_all_alloMaize4Low_16/prunedBy1000/PCA/whole_genome.cov",
                           header = F, stringsAsFactors = F)

pca_1000 <- eigen(cov_data_1000) # take PCA of covariance matrix
n = 4 # make small dataframe w/ only first n eigenvectors
pca_small_1000 <- data.frame(pca_1000$vectors[ , 1:n])
colnames(pca_small_1000) = paste0("PC", 1:n)

# join bams and first 10 PCs of covariance PCA data by position (CAUTION - bam list order and admix results MUST MATCH!)
d_1000 <- bind_cols(pass1_allo4Low, pca_small_1000)  %>%
  #arrange(., popN) %>%
  #arrange(., zea) %>%
  #arrange(., symp_allo) %>%
  mutate(., group = paste(symp_allo, zea, sep = "_"))

# rounded eigen values
PC_var_explained_1000 = round(pca_1000$values, 2)

# plot first PC's
# PC1 and 2, all no filter
p12_coverage_1000 = d_1000 %>%
  ggplot(., aes(PC1, PC2, alpha = .5)) + 
  xlab(paste0("PC1 (", PC_var_explained_1000[1], "%)")) +
  ylab(paste0("PC1 (", PC_var_explained_1000[2], "%)")) +
  ggtitle("PC1 and PC2 (all ind's); every 1000th SNP") +
  geom_point(aes(color = group, size = est_coverage)) +
  scale_colour_manual(values = colors_alphabetical)
ggsave("../plots/PCA_12_by_coverage_pruned1000.png", 
       plot = p12_coverage_1000, device = png(), 
       width = 12, height = 8, units = "in",
       dpi = 200)
# PC1 and 2, excluding super low coverae ind's
p12_1000 = d_1000 %>%
  filter(., est_coverage >= .05) %>%
  ggplot(., aes(PC1, PC2)) + 
  ggtitle("PC1 and PC2 (excl. low coverage ind's < .05x; every 1000th SNP)") +
  xlab(paste0("PC1 (", PC_var_explained_1000[1], "%)")) +
  ylab(paste0("PC1 (", PC_var_explained_1000[2], "%)")) +
  geom_point(aes(color = group)) +
  scale_colour_manual(values = colors_alphabetical)
plot(p12_1000)
ggsave("../plots/PCA_12_pruned1000.png", plot = p12_1000, device = png(), 
       width = 12, height = 8, units = "in",
       dpi = 200)
# PC3 and PC4, excluding low coverage ind's
p34_1000 = d_1000 %>%
  filter(., est_coverage >= .05) %>%
  ggplot(., aes(PC3, PC4)) + 
  ggtitle("PC3 and PC4 (excl. low coverage ind's < .05x); every 1000th SNP") +
  xlab(paste0("PC3 (", PC_var_explained_1000[3], "%)")) +
  ylab(paste0("PC4 (", PC_var_explained_1000[4], "%)")) +
  geom_point(aes(color = group)) +
  scale_colour_manual(values = colors_alphabetical)
plot(p34_1000)
ggsave("../plots/PCA_34_pruned1000.png", plot = p34_1000, device = png(), 
       width = 12, height = 8, units = "in",
       dpi = 200)

# Use vars() to supply faceting variables:
p12_1000 + facet_wrap(~LOCALITY)
ggsave("../plots/PCA_facet_by_population_pruned1000.png", 
       plot = p12_1000 + facet_wrap(~LOCALITY), device = png(), 
       width = 12, height = 8, units = "in",
       dpi = 200)




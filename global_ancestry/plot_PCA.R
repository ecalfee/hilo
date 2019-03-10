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

PREFIX <- "hilo_alloMAIZE_MAIZE4LOW_RIMMA0625_small"
IDs <- data.frame(ID = read.table(paste0("../samples/", PREFIX, "_IDs.list"), header = F, stringsAsFactors = F)$V1, stringsAsFactors = F)
metrics <- read.table(paste0("../filtered_bams/metrics/", PREFIX, ".flagstat.total"), header = F, stringsAsFactors = F)
colnames(metrics) <- c("ID", "total_reads_pass")
hilo <- read.table("../samples/hilo_meta.txt", stringsAsFactors = F, header = T, sep = "\t")
maize4low <- read.table("../samples/MAIZE4LOW_meta.txt", stringsAsFactors = F, header = T)
landraces <- read.table("../samples/alloMAIZE_meta.txt", stringsAsFactors = F, header = T)
seq_Jan2019 <- read.table("../data/HILO_raw_reads/Jan2019_IDs.list", stringsAsFactors = F, header = F)$V1

# a rough coverage estimate is number of reads passing filtering * 150bp/read
# divided by total area they could map to in maize reference genome v4

# size of reference genome reads are mapped to
ref_genome_size <- sum(read.table("../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.fai", 
                                  stringsAsFactors = F)$V2) # V2 is the size of each chromosome mapped to,
# including parts I won't analyze on the Pt and Mt and scaffolds not assigned to chromosomes (but reads would pass Q filters there too)
metrics$est_coverage = round(metrics$total_reads_pass*150/ref_genome_size, 4)

# combine sample meta data
meta <- bind_rows(hilo, maize4low, landraces) %>%
  left_join(., metrics, by = "ID") %>%
  mutate(., group = paste(symp_allo, zea, sep = "_"))
meta$est_coverage[meta$ID=="RIMMA0625"] <- metrics$est_coverage[metrics$ID=="RIMMA0625_small"]
sum(meta$est_coverage < .025)
meta %>%
  filter(., est_coverage < 0.25) %>%
  ggplot(aes(color = group)) +
  geom_boxplot(stat = count)
sum(meta$est_coverage)

# get PCA data
cov_dir <- file.path("results/PCA", PREFIX)
cov_data <- read.table(file.path(cov_dir, "whole_genome.cov"),
                       header = F, stringsAsFactors = F)
# PC's are each column of dataframe pca:
# i.e. PC1 is V1, PC2 is V2 etc.
pca <- eigen(cov_data) # take PCA of covariance matrix
n = 4 # make small dataframe w/ only first n eigenvectors
pca_small <- data.frame(pca$vectors[ , 1:n])
colnames(pca_small) = paste0("PC", 1:n)
# rounded eigen values
PC_var_explained = round(pca$values, 2)
  
# quick plot of new data:
d <- bind_cols(IDs, pca_small) %>%
  left_join(., meta, by = "ID")
d %>%
  filter(., ID %in% seq_Jan2019 | (symp_allo == "allopatric" & zea == "maize")) %>%
  ggplot(., aes(x = -PC1, y = PC2, color = group, 
                size = est_coverage)) + 
  geom_point(alpha = .5) + 
  #scale_colour_manual(values = colors_alphabetical) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC1 (", PC_var_explained[2], "%)")) +
  ggtitle("newly seq. HILO Jan19 only; SNPs spaced >.01cM")

ggsave("plots/PC_new_HILO_seq_Jan19_only.png", 
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)
d %>%
  filter(., !(ID %in% seq_Jan2019)) %>%
  ggplot(., aes(x = -PC1, y = PC2, color = paste(zea, symp_allo, sep = "_"), 
                size = est_coverage)) + 
  geom_point(alpha = .75) + 
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC1 (", PC_var_explained[2], "%)")) +
  ggtitle("old seq. HILO only; SNPs spaced >.01cM")
ggsave("plots/PC_old_HILO_seq_before_Jan19_only.png", 
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)

d %>%
  #filter(., (ID %in% seq_Jan2019 & n <= 200) | (symp_allo == "allopatric" & zea == "maize")) %>%
  filter(., (ID %in% seq_Jan2019 & n <= 200)) %>%
  ggplot(., aes(x = -PC1, y = PC2, color = paste(zea, symp_allo, sep = "_"), 
                size = est_coverage)) + 
  geom_point(alpha = .75) + 
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC1 (", PC_var_explained[2], "%)")) +
  ggtitle("overlap old-new seq. HILO; SNPs spaced >.01cM")


# get metadata for individuals included in NGSadmix analysis
pass1_allo4Low <- read.table("../data/pass1_allo4Low_ids.txt", stringsAsFactors = F, 
                             header = T, sep = "\t")
# join bams and first 10 PCs of covariance PCA data by position (CAUTION - bam list order and admix results MUST MATCH!)
d <- bind_cols(pass1_allo4Low, pca_small)  %>%
  #arrange(., popN) %>%
  #arrange(., zea) %>%
  #arrange(., symp_allo) %>%
  mutate(., group = paste(symp_allo, zea, sep = "_"))

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




#!/usr/bin/env Rscript
# working directory is hilo/
# gets a sample list for individuals with homozygous maize/parv or mexicana alleles
# based on PCA results for putative inversion at mhl1
library(dplyr)
library(tidyr)

# load variables from snakemake
# get covariance matrix estimated by PCAngsd
cov_file = snakemake@input[["cov"]]
# cov_file = "mhl1_inv/results/PCA/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/mhl1_inv.cov"
bam_list = snakemake@input[["bam_list"]]
# bam_list = "samples/HILO_MAIZE55_PARV50_bams.list"

# get output filenames
maize2_parv = snakemake@output[["maize2_parv"]]
#maize2_parv = "mhl1_inv/results/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/maize_cluster_parv_bams.list"
maize2_maize = snakemake@output[["maize2_maize"]]
#maize2_maize = "mhl1_inv/results/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/maize_cluster_maize_bams.list"
mexicana2_maize = snakemake@output[["mexicana2_maize"]]
#mexicana2_maize = "mhl1_inv/results/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/mexicana_cluster_maize_bams.list"
mexicana2_mexicana = snakemake@output[["mexicana2_mexicana"]]
#mexicana2_mexicana = "mhl1_inv/results/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/mexicana_cluster_mexicana_bams.list"

# get sample metadata, including estimated sequence coverage
load(snakemake@input[["meta"]]) 
# load("samples/HILO_MAIZE55_PARV50_meta.RData")
# note: metadata is in the same order as samples input to PCAngsd

# get PCA data
cov_data <- read.table(cov_file, header = F, stringsAsFactors = F)
# PC's are each column of dataframe pca:
# i.e. PC1 is V1, PC2 is V2 etc.
pca <- eigen(cov_data) # take PCA of covariance matrix
m = 2 # make small dataframe w/ only first 2 eigenvectors
# eigenvectors
pca_small <- data.frame(pca$vectors[ , 1:2])
colnames(pca_small) = paste0("PC", 1:2)

# get bams
bams = read.table(bam_list) %>%
  dplyr::rename(bam = V1)

# join ids, bams and first PCs of covariance PCA data by position (CAUTION - bam list order and admix results MUST MATCH!)
d <- bind_cols(meta, bams, pca_small)  %>%
  arrange(., popN) %>%
  arrange(., zea) %>%
  arrange(., symp_allo) %>%
  arrange(., group, ELEVATION)

# cutoffs for maize/parv vs mex clusters based on visual inspection of PC1
maize2_cutoff = 0.02
mexicana2_cutoff = -0.04

# < 5 mexicana in the maize cluster and 
# only 1 parviglumis in the mexicana cluster
# so I won't include these two groups
#d %>% 
#  filter(PC1 > maize2_cutoff) %>%
#  group_by(zea) %>%
#  summarise(n = n())
#d %>% 
#  filter(PC1 < mexicana2_cutoff) %>%
#  group_by(zea) %>%
#  summarise(n = n())

# write files for bams with high confidence assignment to 
# either the maize homozygous allele cluster
# or mexicana homozygous allele cluster
# at the inversion, based on PC1

# parviglumis in the maize cluster
d %>% 
  filter(PC1 > maize2_cutoff & zea == "parviglumis") %>%
  dplyr::select(bam) %>%
  write.table(x = ., 
              file = maize2_parv,
              col.names = F, row.names = F, 
              sep = "\t", quote = F)

# maize in the maize cluster
d %>% 
  filter(PC1 > maize2_cutoff & zea == "maize") %>%
  dplyr::select(bam) %>%
  write.table(x = ., 
              file = maize2_maize,
              col.names = F, row.names = F, 
              sep = "\t", quote = F)

# maize in the mexicana cluster
d %>% 
  filter(PC1 < mexicana2_cutoff & zea == "maize") %>%
  dplyr::select(bam) %>%
  write.table(x = ., 
              file = mexicana2_maize,
              col.names = F, row.names = F, 
              sep = "\t", quote = F)

# mexicana in the mexicana cluster
d %>% 
  filter(PC1 < mexicana2_cutoff & zea == "mexicana") %>%
  dplyr::select(bam) %>%
  write.table(x = ., 
              file = mexicana2_mexicana,
              col.names = F, row.names = F, 
              sep = "\t", quote = F)
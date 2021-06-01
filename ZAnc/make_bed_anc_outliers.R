#!/usr/bin/env Rscript
library(dplyr)

# this script creates a bed file for outlier regions for 
# slope from lm ancestry ~ elev across the genome

# load variables from Snakefile
# zea = "maize"
fdr_file = snakemake@input[["fdr"]]
# fdr_file = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, ".meanAnc.fdr.RData")
anc_file = snakemake@input[["anc"]]
# anc_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/K2/Ne10000_yesBoot/anc/", zea, ".mexicana_anc.bed")
sim_file = snakemake@input[["sim"]]
# sim_file = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, ".MVN.RData")
ancestry = snakemake@params[["ANCESTRY"]]
# ancestry = "mexicana"

# outliers based on FDR:
pos_fdr = snakemake@output[["pos_fdr"]] # FDR high ancestry outliers
# pos_fdr = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, "_pos_mexicana_anc_outliers.fdr05.bed")
neg_fdr = snakemake@output[["neg_fdr"]] # low ancestry outliers
# neg_fdr = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, "_neg_mexicana_anc_outliers.fdr05.bed")

# less stringent than 5% FDR:
# use 5% empirical of the genome (not ancestry tracts) as an outlier cutoff
pos_perc = snakemake@output[["pos_perc"]] # high ancestry outliers
# pos_perc = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, "_pos_mexicana_anc_outliers.perc05.bed")
neg_perc = snakemake@output[["neg_perc"]] # low ancestry outliers
# neg_perc = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, "_neg_mexicana_anc_outliers.perc05.bed")


# use p-value cutoff using from simulated data (p = 0.05 is top 5% simulated points)
pos_p = snakemake@output[["pos_p"]] # high ancestry outliers
# pos_p = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, "_pos_mexicana_anc_outliers.p05.bed")
neg_p = snakemake@output[["neg_p"]] # low ancestry outliers
# neg_p = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, "_neg_mexicana_anc_outliers.p05.bed")

# load data
load(fdr_file) # false discovery rates
load(sim_file) # loads MVN_mean vector
anc <- read.table(anc_file, header = T, stringsAsFactors = F, sep = "\t") # mean ancestry bed file

# write bed files (just outlier tract positions, no header or ancestry)
# filter for positive slope outliers above 5% FDR
anc %>%
  filter(anc_freq > filter(FDRs[[ancestry]], FDR == 0.05 & tail == "high")$threshold) %>%
  dplyr::select(., chr, start, end) %>%
  write.table(., file = pos_fdr,  sep = "\t", quote = F,
              col.names = F, row.names = F) # write bed output file

# filter for low mexicana ancestry outliers
anc %>%
  filter(anc_freq < filter(FDRs[[ancestry]], FDR == 0.05 & tail == "low")$threshold) %>%
  dplyr::select(., chr, start, end) %>%
  write.table(., file = neg_fdr,  sep = "\t", quote = F,
              col.names = F, row.names = F) # write bed output file

# what % of the genome meets the FDR threshold?
print("what % of the genome meets the FDR threshold?")
anc %>%
  filter(anc_freq > filter(FDRs[[ancestry]], FDR == 0.05 & tail == "high")$threshold) %>%
  mutate(length = end - start) %>%
  summarise(total = sum(length)/(2.3*10^9)) %>%
  print(.)

# 5% of genome outliers (less stringent outlier criteria):
# high
anc %>%
  arrange(., anc_freq) %>%
  mutate(length = end - start,
         cum_length = cumsum(length),
         quantile = cum_length/sum(length)) %>%
  filter(quantile > .95) %>%
  arrange(., chr, start) %>%
  dplyr::select(., chr, start, end) %>%
  write.table(., file = pos_perc,  sep = "\t", quote = F,
              col.names = F, row.names = F)
# low
anc %>%
  arrange(., anc_freq) %>%
  mutate(length = end - start,
         cum_length = cumsum(length),
         quantile = cum_length/sum(length)) %>%
  filter(quantile < .05) %>%
  arrange(., chr, start) %>%
  dplyr::select(., chr, start, end) %>%
  write.table(., file = neg_perc,  sep = "\t", quote = F,
              col.names = F, row.names = F)

# p-value from simulations -- p-value = 0.05 cutoff is at quantiles c(2.5%, 97.5%) of simulated points
# high
anc %>%
  filter(anc_freq > quantile(MVN_mean[[ancestry]], 0.975)) %>%
  dplyr::select(., chr, start, end) %>%
  write.table(., file = pos_p,  sep = "\t", quote = F,
              col.names = F, row.names = F)
# low
anc %>%
  filter(anc_freq < quantile(MVN_mean[[ancestry]], 0.025)) %>%
  dplyr::select(., chr, start, end) %>%
  write.table(., file = neg_p,  sep = "\t", quote = F,
              col.names = F, row.names = F)

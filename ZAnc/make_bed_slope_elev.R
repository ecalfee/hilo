#!/usr/bin/env Rscript
library(dplyr)

# this script creates a bed file for outlier regions for 
# slope from lm ancestry ~ elev across the genome

# load variables from Snakefile
# zea = "maize"
fdr_file = snakemake@input[["fdr"]]
# fdr_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".lmElev.fdr.RData")
fit_file = snakemake@input[["fit"]]
# fit_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".lmElev.fit.RData")
sim_file = snakemake@input[["sim"]]
# sim_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".lmElev.sim.RData")
sites_file = snakemake@input[["bed"]]
# sites_file = "local_ancestry/results/thinnedSNPs/HILO_MAIZE55/whole_genome.bed"
fdr_pos = snakemake@output[["fdr_pos"]] # positive slope outliers
# fdr_pos = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, "_pos_lmElev_outliers.fdr05.bed")
fdr_neg = snakemake@output[["fdr_neg"]] # negative slope outliers
# fdr_neg = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, "_neg_lmElev_outliers.fdr05.bed")

# less stringent than 5% FDR:
# use 2% empirical of the genome (not ancestry tracts) as an outlier cutoff
pos_perc = snakemake@output[["pos_perc"]] 
# pos_perc = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, "_pos_lmElev_outliers.perc02.bed")
neg_perc = snakemake@output[["neg_perc"]] 
# neg_perc = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, "_neg_lmElev_outliers.perc02.bed")

# use p-value cutoff using from simulated data (p = 0.05 is top 5% simulated points)
pos_p = snakemake@output[["pos_p"]]
# pos_p = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, "_pos_lmElev_outliers.p05.bed")
neg_p = snakemake@output[["neg_p"]]
# neg_p = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, "_neg_lmElev_outliers.p05.bed")


# load data
load(fdr_file) # false discovery rates
load(fit_file) # observed data: fitted slope mexicana ancestry ~ elevation
load(sim_file) # simulated MVN: fitted slope mexicana ancestry ~ elevation

# site/position information for ancestry tracts
sites_bed <- read.table(sites_file, header = F, stringsAsFactors = F,
                    sep = "\t") %>%
  data.table::setnames(c("chr", "start", "end", "pos"))

# write bed files (just outlier tract positions, no header or slope)
# filter for positive slope outliers above 5% FDR
bind_cols(sites_bed, fits) %>%
  rename(slope = envWeights) %>%
  filter(slope > filter(FDRs, FDR == 0.1 & tail == "high")$thesholds) %>%
  #dplyr::select(., colnames(sites_bed), slope) %>%
  dplyr::select(., chr, start, end) %>%
  write.table(., file = fdr_pos,  sep = "\t", quote = F,
              col.names = F, row.names = F) # write bed output file

# filter for negative slope outliers above 5% FDR
bind_cols(sites_bed, fits) %>%
  rename(slope = envWeights) %>%
  filter(slope < filter(FDRs, FDR == 0.05 & tail == "low")$thesholds) %>%
  dplyr::select(., chr, start, end) %>%
  write.table(., file = fdr_neg,  sep = "\t", quote = F,
              col.names = F, row.names = F) # write bed output file

# what percent of the genome is in an outlier region?
#bind_cols(sites_bed, fits) %>%
#  rename(slope = envWeights) %>%
#  filter(slope > filter(FDRs, FDR == 0.05 & tail == "high")$thesholds) %>%
#  mutate(length = end - start) %>%
#  summarise(total = sum(length)/(2.3*10^9))

# 2% of genome outliers (less stringent outlier criteria):
# high
bind_cols(sites_bed, fits) %>%
  arrange(., envWeights) %>%
  mutate(length = end - start,
         cum_length = cumsum(length),
         quantile = cum_length/sum(length)) %>%
  filter(quantile > .98) %>%
  arrange(., chr, start) %>%
  dplyr::select(., chr, start, end) %>%
  write.table(., file = perc_pos,  sep = "\t", quote = F,
              col.names = F, row.names = F)
# low
bind_cols(sites_bed, fits) %>%
  arrange(., envWeights) %>%
  mutate(length = end - start,
         cum_length = cumsum(length),
         quantile = cum_length/sum(length)) %>%
  filter(quantile < .02) %>%
  arrange(., chr, start) %>%
  dplyr::select(., chr, start, end) %>%
  write.table(., file = perc_neg,  sep = "\t", quote = F,
              col.names = F, row.names = F)

# p-value from simulations -- p-value = 0.05 cutoff is at quantiles c(2.5%, 97.5%) of simulated points
# high
bind_cols(sites_bed, fits) %>%
  filter(envWeights > quantile(fits_sim$envWeights, 0.975)) %>%
  dplyr::select(., chr, start, end) %>%
  write.table(., file = p_pos,  sep = "\t", quote = F,
              col.names = F, row.names = F)
# low
bind_cols(sites_bed, fits) %>%
  filter(envWeights < quantile(fits_sim$envWeights, 0.025)) %>%
  dplyr::select(., chr, start, end) %>%
  write.table(., file = p_neg,  sep = "\t", quote = F,
              col.names = F, row.names = F)
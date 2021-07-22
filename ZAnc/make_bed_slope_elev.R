#!/usr/bin/env Rscript
library(dplyr)

# this script creates a bed file for outlier regions for 
# slope from lm ancestry ~ elev across the genome

# load variables from Snakefile
# zea = "maize"
fdr_file = snakemake@input[["fdr"]]
# fdr_file = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, ".lmElev.fdr.RData")
fit_file = snakemake@input[["fit"]]
# fit_file = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, ".lmElev.fit.RData")
sim_file = snakemake@input[["sim"]]
# sim_file = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, ".lmElev.sim.RData")
sites_file = snakemake@input[["bed"]]
# sites_file = "local_ancestry/results/thinnedSNPs/HILO_MAIZE55/K2/whole_genome.bed"
pos_fdr = snakemake@output[["pos_fdr"]] # positive slope outliers
# pos_fdr = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, "_pos_lmElev_outliers.fdr05.bed")
neg_fdr = snakemake@output[["neg_fdr"]] # negative slope outliers
# neg_fdr = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, "_neg_lmElev_outliers.fdr05.bed")

# less stringent than 5% FDR:
# use 5% empirical of the genome (not ancestry tracts) as an outlier cutoff
pos_perc = snakemake@output[["pos_perc"]] 
# pos_perc = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, "_pos_lmElev_outliers.perc05.bed")
neg_perc = snakemake@output[["neg_perc"]] 
# neg_perc = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, "_neg_lmElev_outliers.perc05.bed")

# use p-value cutoff using from simulated data (p = 0.05 is top 5% simulated points)
pos_p = snakemake@output[["pos_p"]]
# pos_p = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, "_pos_lmElev_outliers.p05.bed")
neg_p = snakemake@output[["neg_p"]]
# neg_p = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, "_neg_lmElev_outliers.p05.bed")


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
  filter(slope > filter(FDRs, FDR == 0.1 & tail == "high")$threshold) %>%
  #dplyr::select(., colnames(sites_bed), slope) %>%
  dplyr::select(., chr, start, end, slope) %>%
  write.table(., file = pos_fdr,  sep = "\t", quote = F,
              col.names = F, row.names = F) # write bed output file

# filter for negative slope outliers above 5% FDR
bind_cols(sites_bed, fits) %>%
  rename(slope = envWeights) %>%
  filter(slope < filter(FDRs, FDR == 0.05 & tail == "low")$threshold) %>%
  dplyr::select(., chr, start, end, slope) %>%
  write.table(., file = neg_fdr,  sep = "\t", quote = F,
              col.names = F, row.names = F) # write bed output file

# what percent of the genome is in an outlier region?
print("what percent of the genome is in an outlier region?")
bind_cols(sites_bed, fits) %>%
  rename(slope = envWeights) %>%
  filter(slope > filter(FDRs, FDR == 0.05 & tail == "high")$threshold) %>%
  mutate(length = end - start) %>%
  summarise(total = sum(length)/(2.3*10^9)) %>%
  print(.)

# 5% of genome outliers (less stringent outlier criteria):
# high
bind_cols(sites_bed, fits) %>%
  arrange(., envWeights) %>%
  mutate(length = end - start,
         cum_length = cumsum(length),
         quantile = cum_length/sum(length)) %>%
  filter(quantile > .95) %>%
  arrange(., chr, start) %>%
  dplyr::select(., chr, start, end, slope) %>%
  write.table(., file = pos_perc,  sep = "\t", quote = F,
              col.names = F, row.names = F)
# low
bind_cols(sites_bed, fits) %>%
  arrange(., envWeights) %>%
  mutate(length = end - start,
         cum_length = cumsum(length),
         quantile = cum_length/sum(length)) %>%
  filter(quantile < .05) %>%
  arrange(., chr, start) %>%
  dplyr::select(., chr, start, end, slope) %>%
  write.table(., file = neg_perc,  sep = "\t", quote = F,
              col.names = F, row.names = F)

# p-value from simulations -- p-value = 0.05 cutoff is at quantiles c(2.5%, 97.5%) of simulated points
# high
bind_cols(sites_bed, fits) %>%
  filter(envWeights > quantile(fits_sim$envWeights, 0.975)) %>%
  dplyr::select(., chr, start, end, slope) %>%
  write.table(., file = pos_p,  sep = "\t", quote = F,
              col.names = F, row.names = F)
# low
bind_cols(sites_bed, fits) %>%
  filter(envWeights < quantile(fits_sim$envWeights, 0.025)) %>%
  dplyr::select(., chr, start, end, slope) %>%
  write.table(., file = neg_p,  sep = "\t", quote = F,
              col.names = F, row.names = F)
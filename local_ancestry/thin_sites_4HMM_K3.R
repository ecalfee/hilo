#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(data.table)

# load variables from Snakefile
source(snakemake@input[["equations"]])
regions_file = snakemake@input[["regions"]]
prefix = snakemake@params[["prefix"]]
min_cM = as.numeric(snakemake@params[["min_cM"]])
min_n_maize = as.numeric(snakemake@params[["min_n_maize"]])
min_n_mex = as.numeric(snakemake@params[["min_n_mex"]])
min_n_parv = as.numeric(snakemake@params[["min_n_parv"]])
top_pbs = as.numeric(snakemake@params[["top_pbs"]])
rdiff_out = snakemake@output[["rdiff"]]
sites_out = snakemake@output[["sites"]]
rpos_out = snakemake@output[["rpos"]]
counts_out = snakemake@output[["counts"]]

# # to test:
# setwd("~/Documents/gitErin/hilo")
#print(getwd())
# source("local_ancestry/FST_PBS_equations.R")
#prefix = "HILO_MAIZE55_PARV50"
#min_cM = 0.001
#min_n_maize = 44
#min_n_mex = 12
#min_n_parv = 40
#top_pbs = 0.10
#rpos_out = "test/TEST3_regions.rpos"
#rdiff_out = "test/TEST3_regions.rdiff"
#sites_out = "test/TEST3_regions.var.sites"
#counts_out = "test/TEST3_counts_thinned_AIMs_pbs.txt"
#regions_file = "test/TEST3_regions.list" # regions 33, 82 and 412

regions = read.table(regions_file, header = F,
                     sep = "\t", stringsAsFactors = F) %>%
  data.table::setnames(., c("chr", "start", "end", "region_n", "txt_file"))

# for each region, takes in a last chr and cM value, and the chromosome and rpos file for a region
# and a minimum cM spacing between thinned sites
# to get a vector of T/F to keep (F = discard) a site
# then uses the vector of T/F to thin the gl file

# read in all the data first
d <- do.call(rbind, lapply(regions$region_n, function(n) {
  rpos0 = read.table(paste0("variant_sites/results/", prefix, "/region_", n, ".rpos"),
                     header = F, sep = "\t", stringsAsFactors = F)$V1
  # load variant sites and reference pop minor allele freqs (maf)
  sites0 = read.table(paste0("variant_sites/results/", prefix, "/region_", n, ".var.sites"),
                      header = F, sep = "\t", stringsAsFactors = F) %>%
    data.table::setnames(c("chr", "pos", "major", "minor"))
  maize_maf = read.table(paste0("local_ancestry/results/alloFreqs/", prefix, "/allopatric_maize/region_", n, ".mafs.gz"),
                         header = T, sep = "\t", stringsAsFactors = F) %>%
    left_join(sites0, ., by = c("chr"="chromo", "pos"="position", "major", "minor")) %>%
    dplyr::mutate(knownEM = ifelse(knownEM == 0, knownEM + 10^-10, knownEM)) %>%
    dplyr::rename(maize_freq = knownEM,
                  maize_nInd = nInd) %>%
    dplyr::mutate(maize_nHap = maize_nInd*2) # high coverage, assume 2 alleles sampled per ind
  mex_maf = read.table(paste0("local_ancestry/results/alloFreqs/", prefix, "/allopatric_mexicana/region_", n, ".mafs.gz"),
                       header = T, sep = "\t", stringsAsFactors = F) %>%
    left_join(sites0, ., by = c("chr"="chromo", "pos"="position", "major", "minor")) %>%
    dplyr::mutate(knownEM = ifelse(knownEM == 0, knownEM + 10^-10, knownEM)) %>%
    dplyr::rename(mex_freq = knownEM,
                  mex_nInd = nInd) %>%
    dplyr::mutate(mex_nHap = mex_nInd) # low coverage, assume 1 allele sampled per ind
  parv_maf = read.table(paste0("local_ancestry/results/alloFreqs/", prefix, "/parv/region_", n, ".mafs.gz"),
                       header = T, sep = "\t", stringsAsFactors = F) %>%
    left_join(sites0, ., by = c("chr"="chromo", "pos"="position", "major", "minor")) %>%
    dplyr::mutate(knownEM = ifelse(knownEM == 0, knownEM + 10^-10, knownEM)) %>%
    dplyr::rename(parv_freq = knownEM,
                  parv_nInd = nInd) %>%
    dplyr::mutate(parv_nHap = parv_nInd*2) # high coverage, assume 2 alleles sampled per ind

  all_maf = bind_cols(mex_maf, 
                      parv_maf[ , c("parv_freq", "parv_nInd", "parv_nHap")], 
                      maize_maf[ , c("maize_freq", "maize_nInd", "maize_nHap")]) %>%
    dplyr::mutate(rpos = rpos0) %>% # position on genetic scale
    # First find SNPs that meet minimum n samples with data to be ancestry informative markers (AIMs)
    dplyr::filter(., !is.na(maize_freq) &
                    !is.na(mex_freq) &
                    !is.na(parv_freq) &
                    maize_nInd >= min_n_maize &
                    mex_nInd >= min_n_mex &
                    parv_nInd >= min_n_parv) %>%
    # calculate differentiation between populations at all SNPs
    dplyr::mutate(
      pbs_maize = sapply(1:nrow(.), function(i)
        pop_branch_stat(p1 = .$maize_freq[i], p2 = .$mex_freq[i], p3 = .$parv_freq[i], 
                        n_hap1 = .$maize_nHap[i], n_hap2 = .$mex_nHap[i], n_hap3 = .$parv_nHap[i])),
      pbs_mex = sapply(1:nrow(.), function(i)
        pop_branch_stat(p1 = .$mex_freq[i], p2 = .$maize_freq[i], p3 = .$parv_freq[i], 
                                              n_hap1 = .$mex_nHap[i], n_hap2 = .$maize_nHap[i], n_hap3 = .$parv_nHap[i])),
      pbs_parv = sapply(1:nrow(.), function(i)
        pop_branch_stat(p1 = .$parv_freq[i], p2 = .$maize_freq[i], p3 = .$mex_freq[i], 
                        n_hap1 = .$parv_nHap[i], n_hap2 = .$maize_nHap[i], n_hap3 = .$mex_nHap[i]))
      )
  }
  ))

# Find SNPs that meet threshold differentiation to be ancestry informative markers (AIMs)

# cutoffs for thinning:
cutoff_maize = quantile(d$pbs_maize, 1 - top_pbs)
cutoff_mex = quantile(d$pbs_mex, 1 - top_pbs)
cutoff_parv = quantile(d$pbs_parv, 1 - top_pbs)

# find AIMs that meet at least one cutoff for thinning
aims = filter(d, pbs_maize >= cutoff_maize | pbs_mex >= cutoff_mex | pbs_parv >= cutoff_parv)

# counts
tot_min_ind = nrow(d)
tot_is_aim = nrow(aims)

# tot_is_aim/tot_min_ind*100 # of SNPs with minimum coverage, percentage that are ancestry informative AIMs

# don't need to keep whole dataframe in memory
rm(d)

# starting values
last_chr = 0
last_cM = -Inf

# which AIM positions to keep?
keep = rep(F, nrow(aims))
for (i in 1:nrow(aims)){
  if (last_chr != aims$chr[i] | aims$rpos[i] - last_cM >= min_cM){
    # keep site if far enough apart (or on a new chromosome)
    keep[i] <- T
    # update last chromosome and cM position for a kept site
    last_chr <- aims$chr[i]
    last_cM <- aims$rpos[i]
  }
}

# SNPs to keep
aims_keep <- aims[keep, ] %>%
  # calculate difference in rpos position in Morgans (divide cM by 100)
  # and set first rdiff of any new chromosome to 1 (only case where current chr will not equal previous chr)
  dplyr::mutate(rdiff = ifelse(chr == lag(chr, default = 0), (rpos - lag(rpos, default = 0))/100, 1)) 

#aims_keep %>%
#  filter(rdiff != 1) %>%
#  summarise(median = median(rdiff)*100,
#            mean = mean(rdiff)*100,
#            min = min(rdiff)*100,
#            max = max(rdiff)*100)

# write output files (append to file if it's not the first region)
options(scipen = 999) # do not print scientific notation
write.table(aims_keep$rpos, file = rpos_out, append = F, col.names = F, row.names = F, quote = F, sep = "\t")
write.table(dplyr::select(aims_keep, chr, pos, major, minor), file = sites_out, append = F, col.names = F, row.names = F, quote = F, sep = "\t")
write.table(aims_keep$rdiff, file = rdiff_out, append = F, col.names = F, row.names = F, quote = F, sep = "\t")

# write summary of total counts:
tot_keep = sum(keep)
data.frame(variable = c("min_ind", "is_aim", "keep", "top_pbs", "min_pbs_maize", "min_pbs_mex", "min_pbs_parv", "perc_pbs_maize", "perc_pbs_mex", "perc_pbs_parv"),
           value = c(tot_min_ind, tot_is_aim, tot_keep, top_pbs, cutoff_maize, cutoff_mex, cutoff_parv, 
                      sum(aims_keep$pbs_maize >= cutoff_maize)/nrow(aims_keep)*100, sum(aims_keep$pbs_mex >= cutoff_mex)/nrow(aims_keep)*100, sum(aims_keep$pbs_parv >= cutoff_parv)/nrow(aims_keep)*100),
           stringsAsFactors = F) %>%
  format(., digits = 2, scientific = F) %>%
  write.table(., file = counts_out, col.names = T, row.names = F, quote = F, sep = "\t")
options(scipen = 0) # back to default

#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(data.table)


# load variables from Snakefile
regions_file = snakemake@input[["regions"]]
prefix_all = snakemake@params[["prefix_all"]]
min_cM = as.numeric(snakemake@params[["min_cM"]])
min_n_maize = as.numeric(snakemake@params[["min_n_maize"]])
min_n_mex = as.numeric(snakemake@params[["min_n_mex"]])
min_maf_diff = as.numeric(snakemake@params[["min_maf_diff"]])
rdiff_out = snakemake@output[["rdiff"]]
sites_out = snakemake@output[["sites"]]
rpos_out = snakemake@output[["rpos"]]
counts_out = snakemake@output[["counts"]]

# # to test:
<<<<<<< HEAD
#setwd("~/Documents/gitErin/hilo")
=======
# setwd("~/Documents/gitErin/hilo")
>>>>>>> 1ccec06eef065c379a427dde73dd638e7d7b25c8
#print(getwd())
#prefix_all = "HILO_MAIZE55"
#min_cM = 0.001
#min_n_maize = 44
#min_n_mex = 12
#min_maf_diff = 0.3
#rpos_out = "test/TEST2_regions.rpos"
#rdiff_out = "test/TEST2_regions.rdiff"
#sites_out = "test/TEST2_regions.var.sites"
#counts_out = "test/counts_thinned_AIMs.txt"
#regions_file = "test/TEST2_regions.list"

regions = read.table(regions_file, header = F,
                     sep = "\t", stringsAsFactors = F) %>%
  data.table::setnames(., c("chr", "start", "end", "region_n", "txt_file"))

<<<<<<< HEAD
=======
#print(regions)
>>>>>>> 1ccec06eef065c379a427dde73dd638e7d7b25c8

# for each region, takes in a last chr and cM value, and the chromosome and rpos file for a region
# and a minimum cM spacing between thinned sites
# to get a vector of T/F to keep (F = discard) a site
# then uses the vector of T/F to thin the gl file

# keep count of sites meeting inclusion thresholds
tot_min_ind = 0
tot_min_diff = 0
tot_is_aim = 0       
tot_keep = 0

# starting values
last_chr = 0
last_cM = -Inf

for (j in 1:nrow(regions)){
  chr = regions$chr[j] # each region spans only 1 chromosome
  n = regions$region_n[j] # which region

  rpos0 = read.table(paste0("variant_sites/results/", prefix_all, "/region_", n, ".rpos"), 
                    header = F, sep = "\t", stringsAsFactors = F)$V1
  # load variant sites and reference pop minor allele freqs (maf)
  sites0 = read.table(paste0("variant_sites/results/", prefix_all, "/region_", n, ".var.sites"),
                     header = F, sep = "\t", stringsAsFactors = F) %>%
    data.table::setnames(c("chr", "pos", "major", "minor"))
  maize_maf = read.table(paste0("variant_sites/results/popFreq/allopatric_maize/region_", n, ".mafs.gz"), 
                         header = T, sep = "\t", stringsAsFactors = F) %>%
    left_join(sites0, ., by = c("chr"="chromo", "pos"="position", "major", "minor"))
  mex_maf = read.table(paste0("variant_sites/results/popFreq/allopatric_mexicana/region_", n, ".mafs.gz"), 
                         header = T, sep = "\t", stringsAsFactors = F) %>%
    left_join(sites0, ., by = c("chr"="chromo", "pos"="position", "major", "minor"))
 
  # First find SNPs that meet threshold difference in allele frequency 
  # and minimum n samples with data to be ancestry informative markers (AIMs)
  min_ind = !is.na(maize_maf$phat) & !is.na(mex_maf$phat) &
               maize_maf$nInd >= min_n_maize & mex_maf$nInd >= min_n_mex
  min_diff = !is.na(abs(maize_maf$phat - mex_maf$phat)) &
    (abs(maize_maf$phat - mex_maf$phat) >= min_maf_diff)
  is_aim = min_ind & min_diff  
  rpos_aims = rpos0[is_aim]
  sites_aims = sites0[is_aim, ]
    
  # second, thin AIMs  
  # which positions to keep?
  keep = rep(F, length(rpos_aims))
  for (i in 1:length(rpos_aims)){
    #print(paste("i:", i, "last_chr:", last_chr, "chr:", chr, "rpos_aims:", rpos_aims[i], "last_cM:", last_cM, "min_cM:", min_cM))
    if (last_chr != chr | rpos_aims[i] - last_cM >= min_cM){
      # keep site if far enough apart (or on a new chromosome)
      keep[i] <- T
      # update last chromosome and cM position for a kept site
      last_chr <- chr
      last_cM <- rpos_aims[i]
    }
  }
  
  # add to counts
  tot_min_ind = tot_min_ind + sum(min_ind)
  tot_min_diff = tot_min_diff + sum(min_diff)
  tot_is_aim = tot_is_aim + sum(is_aim)
  tot_keep = tot_keep + sum(keep)
  
  # SNPs to keep
  sites_keep <- sites_aims[keep, ]
  rpos_keep <- rpos_aims[keep]
  
  # calculate difference in rpos position in Morgans (divide cM by 100)
  rdiff <- c(0, diff(rpos_keep/100))
  # set first rdiff of each chromosome = 1
  rdiff[c(T, diff(sites_keep$chr) != 0)] <- 1
  
  # write output files (append to file if it's not the first region)
  options(scipen = 999) # do not print scientific notation
  write.table(rpos_keep, file = rpos_out, append = (n > 0), col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(sites_keep, file = sites_out, append = (n > 0), col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(rdiff, file = rdiff_out, append = (n > 0), col.names = F, row.names = F, quote = F, sep = "\t")
  options(scipen = 0) # back to default
}

# write summary of total counts:
data.frame(filter = c("min_ind", "min_maf_diff", "is_aim", "all"), 
           n_pass = c(tot_min_ind, tot_min_diff, tot_is_aim, tot_keep),
           stringsAsFactors = F) %>%
write.table(., file = counts_out, col.names = T, row.names = F, quote = F, sep = "\t")


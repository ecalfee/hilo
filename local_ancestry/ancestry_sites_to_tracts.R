#!/usr/bin/env Rscript

# this script converts var sites used in hmm ancestry inference
# into tracts around those sites
# specifially, converts chr and pos to a bed file, where start and end are midway between
# the original snp positions with ancestry calls, and the 'value' of the tract is the original ancestry call snp position
# (not necessarily the center of the tract because of asymmetry in snp spacing)
# first and last ancestry call have tracts that end at that call (rather than extending to the chr end)

# to run:
# Rscript ancestry_sites_to_tracts.R pass2_alloMAIZE

library(dplyr)
library(tidyr)

# genomic sites information
# SNPs with ancestry calls
args = commandArgs(trailingOnly=TRUE)
PREFIX = args[1]
# load sites
dir_sites = paste0("../local_ancestry/results/thinnedSNPs/", PREFIX)
sites <- do.call(rbind, 
                 lapply(1:10, function(i)
                   read.table(paste0(dir_sites, "/chr", i, ".var.sites"),
                              header = F, stringsAsFactors = F)))
colnames(sites) <- c("chr", "pos", "major", "minor")
# convert to tracts
sites2 <- head(sites)
tracts <- sites %>%
  mutate(end = pos + round(c(diff(pos), 0)/2) - 1,
         start = c(0, end[1:(length(end) - 1)]),
         first = c(1, diff(chr)) == 1, # first pos on chr
         last = c(diff(chr), 1) == 1, # last pos on chr
         end = ifelse(last, pos, end),
         start = ifelse(first, pos - 1, start)) %>%
  dplyr::select(chr, start, end, pos) # nope!
write.table(tracts,
            paste0(dir_sites, "/tracts.bed"), 
            quote = F, col.names = F, row.names = F, 
            sep = "\t")

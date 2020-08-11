#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)

# print working directory
print(paste("R working directory:", getwd()))

# load variables from Snakefile
rmap_file = snakemake@input[["rmap"]] # recombination map
# rmap_file = "../data/linkage_map/ogut_fifthcM_map_agpv4_INCLUDE.txt"
sites_file = snakemake@input[["sites"]]
# sites_file = "../variant_sites/results/TEST/region_1.var.sites"
out_file = snakemake@output[["rpos"]]
# out_file = "../variant_sites/results/TEST/region_1.rpos"

# load sites file
sites = read.table(sites_file, header = F, sep = "\t") %>%
  data.table::setnames(., c("chr", "pos", "allele1", "allele2"))

# load recombination map
rmap = read.table(rmap_file, header = F, sep = "\t") %>%
  data.table::setnames(., c("sts", "marker", "pos_cM", "chr", "pos")) %>%
  filter(chr == unique(sites$chr)) # filter to just relevant chromosome

# interpolate/extrapolate map positions (cM) for all sites using a monotone cubic spline
sites$rpos = spline(x = rmap$pos, y = rmap$pos_cM, 
              xout = sites$pos, method = "hyman")$y

# write output file
write.table(sites$rpos, out_file, 
            col.names = F, row.names = F, 
            sep = "\t", quote = F)

# quick plot to check
# plot(rmap$pos, rmap$pos_cM)
# points(sites$pos, sites$rpos, col = "blue")
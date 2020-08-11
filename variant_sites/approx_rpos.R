#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)

# print working directory
# setwd("~/Documents/gitErin/hilo")
print(paste("R working directory:", getwd()))

# load variables from Snakefile
rmap_file = snakemake@input[["rmap"]] # recombination map
# rmap_file = "data/linkage_map/ogut_fifthcM_map_agpv4_EXTENDED.txt"
sites_file = snakemake@input[["sites"]]
# sites_file = "variant_sites/results/TEST/region_1.var.sites"
out_file = snakemake@output[["rpos"]]
# out_file = "variant_sites/results/TEST/region_1.rpos"

# load sites file
sites = read.table(sites_file, header = F, sep = "\t") %>%
  data.table::setnames(., c("chr", "pos", "allele1", "allele2"))

# load recombination map
rmap = read.table(rmap_file, header = T, sep = "\t", stringsAsFactors = F) %>%
  filter(chr == unique(sites$chr))


# interpolate map positions (cM) for all sites using linear approximation
# from Ogut 0.2cM genetic map
# note: for points outside of the original linkage map (at ends of the chr),
# we use the recombination rate from the nearest mapped window ('EXTENDED' map)
sites$rpos = stats::approx(x = rmap$pos_bp, y = rmap$pos_cM,
              xout = sites$pos, method = "linear")$y

# write output file
write.table(sites$rpos, out_file,
            col.names = F, row.names = F,
            sep = "\t", quote = F)

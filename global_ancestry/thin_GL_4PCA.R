#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)

# load variables from Snakefile
prefix_all = snakemake@params[["prefix_all"]]
min_cM = snakemake@params[["min_cM"]]
file_out = snakemake@output[["gl"]]

regions_file = snakemake@input[["regions"]]
# regions_file = "data/refMaize/divide_5Mb/ALL_regions.list"

regions = read.table(regions_file, header = F,
                     sep = "\t", stringsAsFactors = F) %>%
  data.table::setnames(., c("chr", "start", "end", "region_n", "txt_file"))

# for each region, takes in a last chr and cM value, and the chromosome and rpos file for a region
# and a minimum cM spacing between thinned sites
# to get a vector of T/F to keep (F = discard) a site
# then uses the vector of T/F to thin the gl file

# starting values
last_chr = 0
last_cM = NA

for (j in 1:nrow(regions)){
  chr = regions$chr[j] # each region spans only 1 chromosome
  n = regions$region_n[j] # which region
  rpos_file = paste0("variant_sites/results/", prefix_all, "/region_", regions$region_n[j], ".rpos")
  rpos = read.table(rpos_file, 
                    header = F, sep = "\t", stringsAsFactors = F)
  # which positions to keep?
  keep = rep(F, length(rpos))
  for (i in 1:length(rpos)){
    if (last_chr != sites$chr[i] | rpos[i] - last_cM >= min_cM){
      # keep site if far enough apart (or on a new chromosome)
      keep[i] <- T
      # update last chromosome and cM position for a kept site
      last_chr <- sites$chr[i]
      last_cM <- rpos[i]
    }
  }
  # thin gl file:
  gl_file = paste0("variant_sites/results/", prefix_all, "/region_", regions$region_n[j], ".beagle.gz")
  thin_gl = read.table(gzfile(gl_file),
                  check.names = F, # no checknames because there are repeated header values in beagle format
                  header = T, stringsAsFactors = F, 
                  sep = "\t")[keep, ] # keep only a subset of rows
  
  # concat thinned GL's to output file (compressed)
  write.table(thin_gl, gzfile(file_out),
              append = n > 0,
              col.names = n == 0, # only write column names if it's region 0
              quote = F, sep = "\t",
              col.names = T, row.names = F)
}
#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(data.table)

# load variables from Snakefile
prefix_all = snakemake@params[["prefix_all"]]
min_cM = as.numeric(snakemake@params[["min_cM"]])
file_out = snakemake@output[["gl"]]
regions_file = snakemake@input[["regions"]]

# # to test:
# setwd("~/Documents/gitErin/hilo")
# prefix_all = "HILO_MAIZE55"
# min_cM = 0.01
# file_out = "test/whole_genome.beagle.gz"
# regions_file= "test/TEST2_regions.list"

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
  rpos_file = paste0("variant_sites/results/", prefix_all, "/region_", n, ".rpos")
  rpos = read.table(rpos_file, 
                    header = F, sep = "\t", stringsAsFactors = F)$V1
  # which positions to keep?
  keep = rep(F, length(rpos))
  for (i in 1:length(rpos)){
    if (last_chr != chr | rpos[i] - last_cM >= min_cM){
      # keep site if far enough apart (or on a new chromosome)
      keep[i] <- T
      # update last chromosome and cM position for a kept site
      last_chr <- chr
      last_cM <- rpos[i]
    }
  }
  # thin gl file:
  gl_file = paste0("variant_sites/results/", prefix_all, "/region_", n, ".beagle.gz")
  
  thin_gl <- data.table::fread(gl_file)[keep, ] # keep only a subset of rows

  # concat thinned GL's to output file (compressed)
  data.table::fwrite(thin_gl, 
                     file_out,
                     sep = "\t", quote = F, 
                     compress = "gzip",
                     row.names = F,
                     col.names = (n == 0),
                     append = (n > 0))
  
  # # fread is a faster version of read.table and write.table
  # thin_gl = read.table(gzfile(gl_file),
  #                      check.names = F, # no checknames because there are repeated header values in beagle format
  #                      header = T, stringsAsFactors = F, 
  #                      sep = "\t")[keep, ] # keep only a subset of rows
  # 
  # write.table(thin_gl, gzfile(file_out),
  #             append = (n > 0), # if region 0 (first region) start new file, otherwise append
  #             col.names = (n == 0), # only write column names if it's region 0
  #             row.names = F,
  #             quote = F, sep = "\t")
}

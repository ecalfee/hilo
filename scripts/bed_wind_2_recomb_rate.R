#!/usr/bin/env Rscript

# this script takes in a bed file and calculate the local recombination rate for
# each window in that bed file, using linear interpolation

library(dplyr)
library(tidyr)
source("calcMapLinearApprox.R") # loads recomb map and function bp2cM()

# to run:
# Rscript bed_wind_2_recomb_rate.R ../data/refMaize/windows_10kb/whole_genome.bed

# arguments
args = commandArgs(trailingOnly=TRUE)
# input (bed file):
input_file = args[1]
# output is same file name with .recomb instead of .bed extension:
output_file = paste0(tools::file_path_sans_ext(input_file), ".recomb")

# read input.bed
bed = read.table(input_file, sep = "\t", header = F, stringsAsFactors = F)[ , 1:3]
colnames(bed) = c("chr", "start", "end")

# calculate distance in cM and in bp between start and end of window
dist_cM = apply(bed, 1, function(row)
  bp2cM(chrom = row["chr"], pos_bp = row["end"]) - bp2cM(chrom = row["chr"], pos_bp = row["start"] + 1))
dist_bp = bed$end - (bed$start + 1)

r = dist_cM/dist_bp * 10^6 # recomb rate in cM/Mb units

# write output file with recombination rates per window
write.table(r, output_file, col.names = F, row.names = F, quote = T)
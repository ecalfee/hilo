#!/usr/bin/env Rscript

# input:
# chromosome #, hilo individual #, and directory where the countsACGT subdirectory is located
# this script finds ACGT read countsfile for an individual at a chromosome:
# dir/countsACGT/chrI/hilo_X.counts.gz
# positions with counts (>0 coverage) hilo_X.pos.gz
# and it finds the dir/chrI.var.sites file 
# with major/minor allele and all positions desired for counts

# output: file with two columns: read counts for major allele,
# and read counts for minor allele for all positions
# (even for positions where the individual has zero coverage, ie. 0 0 counts)
# each individual for each chromosome has an output file
# the same input directory will be used for output: 
# dir/countsMajMin/chrI/hilo_X.counts.txt

library(dplyr)
# to run:
# Rscript count_reads_major_minor.R 10 12 ../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/

print(getwd()) # print current directory

# arguments
args = commandArgs(trailingOnly=TRUE)
# chromosome #
#CHR = 10
CHR = as.integer(args[1])
# hilo id #
ID = as.integer(args[2])
# path to directory
#path = "../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM"
path = args[3]

acgt_file = paste0(path, "/countsACGT/chr", CHR, "/hilo_", ID, ".counts.gz")
pos_file = paste0(path, "/countsACGT/chr", CHR, "/hilo_", ID, ".pos.gz")
sites_file = paste0(path, "/chr", CHR, ".var.sites")
output_directory = paste0(path, "/countsMajMin/chr", CHR)
output_file = paste0(output_directory, "/hilo_", ID, ".counts.txt")

# input data
SNPs = read.table(sites_file,
                  header = F, stringsAsFactors = F, sep = "\t")
colnames(SNPs) = c("chr", "pos", "major", "minor")
acgt = read.table(acgt_file, stringsAsFactors = F, header = T)
colnames(acgt) = substr(colnames(acgt), 4, 4) # make totA -> A column names
pos = read.table(pos_file, stringsAsFactors = F, header = T)
d = cbind(pos, acgt) %>%
  left_join(SNPs, ., by = c("chr", "pos")) %>%
  tidyr::gather(., "allele", "n", c("A", "C", "G", "T"))
maj_counts = filter(d, major==allele) %>%
  select(-allele) # get read counts for alleles matching major allele
min_counts = filter(d, minor==allele) %>%
  select(-allele) # get read counts for allele matching minor allele
counts = full_join(maj_counts, min_counts, 
                   suffix = c("_major", "_minor"),
                   by = c("chr", "pos", "major", "minor", "totDepth")) %>%
  .[order(.$pos),] # re-order by nucleotide position within chromosome

#with(counts, sum(is.na(n_major) & !is.na(n_minor)))
# make output directory if it doesn't exist
if (!dir.exists(output_directory)) dir.create(output_directory, recursive = T) # will simply warn if it directory already exists

# write output file
write.table(select(counts, c(n_major, n_minor)), 
            output_file,
            na = "0", # write NA's (no coverage) as zero counts
            row.names = F, col.names = T, quote = F) # no headers, just counts
warnings() # print any warnings

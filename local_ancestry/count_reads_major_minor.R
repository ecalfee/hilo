#!/usr/bin/env Rscript

# input:
# chromosome #, hilo individual #, and directory where the countsACGT subdirectory is located
# this script finds ACGT read counts file for an individual at a chromosome:
# DIR_IN/chr1.counts.gz
# positions with counts (>0 coverage) DIR_IN/chr1.pos.gz
# and it finds the DIR_SITES/chr1.var.sites file
# with major/minor allele and all positions desired for counts

# output: file with two columns: read counts for major allele,
# and read counts for minor allele for all positions
# (even for positions where the individual has zero coverage, ie. 0 0 counts)
# each individual for each chromosome has a diff. output file
# DIR_OUT/chr1.counts.txt

library(dplyr)
# to run:
# Rscript count_reads_major_minor.R 10 HILO12 results/thinnedSNPs/pass2_alloMAIZE results/countsACGT/pass2_alloMAIZE/HILO12 results/countsMajMin/pass2_alloMAIZE/HILO12

print(getwd()) # print current directory

# arguments
args = commandArgs(trailingOnly=TRUE)
# chromosome #
#CHR = 10
CHR = as.integer(args[1])
# sample ID
ID = args[2]
# path to directories
DIR_SITES = args[3]
DIR_IN = args[4]
DIR_OUT = args[5]

acgt_file = paste0(DIR_IN, "/chr", CHR, ".counts.gz")
pos_file = paste0(DIR_IN, "/chr", CHR, ".pos.gz")
sites_file = paste0(DIR_SITES, "/chr", CHR, ".var.sites")
output_file = paste0(DIR_OUT, "/chr", CHR, ".counts.txt")

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

# write output file
write.table(select(counts, c(n_major, n_minor)),
            output_file,
            na = "0", # write NA's (no coverage) as zero counts
            row.names = F, col.names = T, quote = F) # no headers, just counts
warnings() # print any warnings

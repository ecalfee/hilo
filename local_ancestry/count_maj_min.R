#!/usr/bin/env Rscript
library(dplyr)

# this script converts ACGT read counts into major/minor allele read counts

# input: ACGT read counts file for an individual *.counts.gz,
# positions with counts (>0 coverage) *.pos.gz,
# and it finds the *.var.sites file with major/minor allele and all positions desired for counts

# output: file with two columns: read counts for major allele,
# and read counts for minor allele for all positions *.counts.txt
# (even for positions where the individual has zero coverage, ie. 0 0 counts)

# load variables from Snakefile
sites_file = snakemake@input[["sites"]]
acgt_file = snakemake@input[["acgt"]]
pos_file = snakemake@input[["pos"]]
output_file = snakemake@output[["majmin"]]

# to test:
#ID = "HILO9"
#prefix_all = "HILO_MAIZE55"
#sites_file = paste0("local_ancestry/results/thinnedSNPs/", prefix_all, "whole_genome.var.sites")
#acgt_file = paste0("local_ancestry/results/countsACGT/", prefix_all, "/", ID, ".counts.gz")
#pos_file = paste0("local_ancestry/results/countsACGT/", prefix_all, "/", ID, ".pos.gz")
#output_file = paste0("local_ancestry/results/countsMajMin/", prefix_all, "/", ID, ".counts.txt")

# get data
SNPs = read.table(sites_file,
                  header = F, stringsAsFactors = F, sep = "\t")
colnames(SNPs) = c("chr", "pos", "major", "minor")
acgt = read.table(acgt_file, stringsAsFactors = F, header = T)
colnames(acgt) = substr(colnames(acgt), 4, 4) # make totA -> A column names
pos = read.table(pos_file, stringsAsFactors = F, header = T)
d = cbind(pos, acgt) %>%
  left_join(SNPs, ., by = c("chr", "pos")) %>%
  tidyr::gather(., "allele", "n", c("A", "C", "G", "T"))

# do counts
maj_counts = filter(d, major==allele) %>%
  select(-allele) # get read counts for alleles matching major allele
min_counts = filter(d, minor==allele) %>%
  select(-allele) # get read counts for allele matching minor allele
counts = full_join(maj_counts, min_counts,
                   suffix = c("_major", "_minor"),
                   by = c("chr", "pos", "major", "minor", "totDepth")) %>%
  arrange(chr, pos) # re-order by nucleotide position within chromosomes

# write output file
write.table(select(counts, c(n_major, n_minor)),
            output_file,
            na = "0", # write NA's (no coverage) as zero counts
            row.names = F, col.names = T, quote = F) # no headers, just counts
warnings() # print any warnings

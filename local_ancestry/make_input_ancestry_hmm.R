#!/usr/bin/env Rscript
library(dplyr)

# this script takes in a hilo sympatric population
# a file produced by make_allo_counts_ancestry_hmm.R with allopatric counts
# and SNP information, and genetic map distance between markers
# and outputs counts and ploidy files ready for ancestry_hmm analysis
# also outputs a simple list, in order, of the admix individuals in the input files

# load variables from Snakefile
prefix_all = snakemake@params[["prefix_all"]]
allo_file = snakemake@input[["allo"]]
pop_ids_file = snakemake@input[["pop_ids"]]
output_file_counts = snakemake@output[["counts"]]
output_file_ploidy = snakemake@output[["ploidy"]]

# to test:
# setwd(~/Documents/gitErin/hilo)
# prefix_all = "HILO_MAIZE55"
# allo_file = paste0("local_ancestry/results/thinnedSNPs/", prefix_all, "/whole_genome.allo.counts")
# POP = "pop360" # just an example
# pop_ids_file = paste0("samples/Over0.5x_byPop/", POP, "_ids.list")
# output_file_counts = paste0("local_ancestry/results/ancestry_hmm/", prefix_all + "/input/", POP, ".counts")
# output_file_ploidy = paste0("local_ancestry/results/ancestry_hmm/", prefix_all + "/input/", POP, ".ploidy")

POP_IDs = read.table(pop_ids_file, stringsAsFactors = F, header = F)$V1
DIR_COUNTS = paste0("local_ancestry/results/countsMajMin/", prefix_all)
ALLO_COLS = read.table(allo_file, header = T, stringsAsFactors = F, sep = " ") # space separated file

# get counts for each sample id in the population
POP_COUNTS = lapply(POP_IDs, function(id) 
  read.table(paste0(DIR_COUNTS, "/", id, ".counts.txt"), 
             header = T, stringsAsFactors = F) %>%
    data.table::setnames(paste0(id, "_", colnames(.)))) # label columns with the sample id and major/minor

options(scipen = 999) # don't write output in scientific notation
# combine data for counts output file
bind_cols(ALLO_COLS, POP_COUNTS) %>%
  write.table(., output_file_counts,
              na = "0", # write NAs as zero (no ref or alt allele counts for that individual)
              row.names = F, col.names = F, quote = F, sep = " ")

# write ploidy file for included individuals in the same order they appear in ancestry_hmm input file
data.frame(id = POP_IDs, ploidy = 2, stringsAsFactors = F) %>%
  write.table(., output_file_ploidy,
            row.names = F, col.names = F, quote = F, sep = " ")



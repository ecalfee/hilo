#!/usr/bin/env Rscript

# this script takes in a hilo admixed population #
# a file produced by make_allo_counts_ancestry_hmm.R with allopatric counts
# and SNP information, and genetic map distance between markers
# and outputs an input file ready for ancestry_hmm analysis
# also outputs a simple list, in order, of the admix individuals in the input file

library(dplyr)

# to run:
# Rscript make_input_ancestry_hmm.R 360 ../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM ../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/input

# arguments
args = commandArgs(trailingOnly=TRUE)
# population number
N = as.integer(args[1])
# paths to input and output directories
dir_input = args[2]
dir_output = args[3]
if (!file.exists(dir_output)){ # make sure output directory exists first (or create)!
  dir.create(file.path(dir_output), recursive = T)
}

min_coverage = 0.05 # minimum coverage for an individual to be included
# in local ancestry assessment

# find allopatric individuals to include
pass1 <- read.table("../data/pass1_ids.txt", stringsAsFactors = F,
                    header = T, sep = "\t")
# find individuals from the current population AND that pass minimum coverage
pop_ids <- filter(pass1, popN == N & est_coverage >= min_coverage)$n

for (i in 1:10){ # for each chromosome, get counts
  # from allopatric reference panels
  allo_counts = read.table(paste0(dir_input, "/allo_counts_chr", i, ".txt"),
                           header = T, stringsAsFactors = F)
  # and each admixed individual
  d = allo_counts
  for (id in pop_ids){
    ind_counts_file = paste0(dir_input, "/countsMajMin/chr", i, "/hilo_", id, ".counts.txt")
    ind_counts_new = read.table(ind_counts_file, header = T, stringsAsFactors = F, sep = "\t")
    colnames(ind_counts_new) = paste0(id, "_", colnames(ind_counts_new))
    d = cbind(d, ind_counts_new) # column bind b/c all files have the same sites
  }
  # write lines for current chromosome to file
  select(d, -major, -minor) %>%
  write.table(., paste0(dir_output, "/pop", N, ".anc_hmm.input"),
              na = "0", # write NAs as zero (no ref or alt allele counts for that individual)
              append = (i != 1), # chr 1 creates a new file, later chromosomes append
              row.names = F, col.names = F, quote = F)
}
# write id's of included individuals in same order they appear in ancestry_hmm input file
write.table(pop_ids, paste0(dir_output, "/pop", N, ".anc_hmm.ids"),
            row.names = F, col.names = F, quote = F)

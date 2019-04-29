#!/usr/bin/env Rscript

# this script takes in a hilo admixed population #
# a file produced by make_allo_counts_ancestry_hmm.R with allopatric counts
# and SNP information, and genetic map distance between markers
# and outputs an input file ready for ancestry_hmm analysis
# also outputs a simple list, in order, of the admix individuals in the input file

library(dplyr)

# to run:
# Rscript make_input_ancestry_hmm.R 360 results/counts/pass2_alloMAIZE results/countsMajMin/pass2_alloMAIZE ../global_ancestry/results/NGSAdmix/pass1_alloMAIZE/globalAdmixtureByIncludedIndividual.txt results/ancestry_hmm/pass2_alloMAIZE/input

# arguments
args = commandArgs(trailingOnly=TRUE)
# population number
N = as.integer(args[1])
# paths to input and output directories
dir_allo_counts = args[2]
dir_symp_counts = args[3]
inds_file = args[4] # has list of all included individuals, and which pop
dir_output = args[5]

# which individuals are included and part of this pop?
hilo <- read.table(inds_file, stringsAsFactors = F,
                    header = T, sep = "\t")
# find individuals from the current population AND that pass minimum coverage
pop_ids <- filter(hilo, popN == N)$ID

for (i in 1:10){ # for each chromosome, get counts
  # from allopatric reference panels
  allo_counts = read.table(paste0(dir_allo_counts, "/allo_counts_chr", i, ".txt"),
                           header = T, stringsAsFactors = F)
  # and each admixed individual
  for (id in pop_ids){
    ind_counts_file = paste0(dir_symp_counts, "/", id, "/chr", i, ".counts.txt")
    ind_counts_new = read.table(ind_counts_file, header = T, stringsAsFactors = F, sep = "\t")
    colnames(ind_counts_new) = paste0(id, "_", colnames(ind_counts_new))
    d = cbind(allo_counts, ind_counts_new) # column bind b/c all files have the same sites
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
# write ploidy file for included individuals in the same order they appear in ancestry_hmm input file
ploidy = data.frame(id = paste0("HILO", pop_ids), ploidy = 2)
write.table(ploidy, paste0(dir_output, "/pop", N, ".anc_hmm.ids.ploidy"),
            row.names = F, col.names = F, quote = F, sep = "\t")



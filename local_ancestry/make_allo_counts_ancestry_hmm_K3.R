#!/usr/bin/env Rscript
library(dplyr)

# output:
# this script makes a file with allopatric allele counts
# and SNP information, and genetic map distance between markers
# this output is used by make_input_ancestry_hmm.R, which adds population-specific
# read counts for the admixed population being run

print(getwd()) # print current directory

# load variables from Snakefile
sites_file = snakemake@input[["sites"]]
rdiff_file = snakemake@input[["rdiff"]]
mex_file = snakemake@input[["mex"]]
parv_file = snakemake@input[["parv"]]
maize_file = snakemake@input[["maize"]]
output_file = snakemake@output[["allo_counts"]]
DIR_COUNTS = snakemake@params[["dir_counts"]]

# to test:
# setwd(~/Documents/gitErin/hilo)
# prefix_all = "HILO_MAIZE55_PARV50/K3"
# sites_file = paste0("local_ancestry/results/thinnedSNPs/", prefix_all, "/whole_genome.var.sites")
# rdiff_file = paste0("local_ancestry/results/thinnedSNPs/", prefix_all, "/whole_genome.rdiff")
# mex_file = "samples/ALL_byPop/allopatric_mexicana_ids.list"
# parv_file = "samples/ALL_byPop/allopatric_parviglumis_ids.list"
# maize_file = "samples/ALL_byPop/allopatric_maize_ids.list"
# output_file = paste0("local_ancestry/results/thinnedSNPs/", prefix_all, "/whole_genome.allo.counts")

# input data
SNPs = read.table(sites_file,
                  header = F, stringsAsFactors = F, sep = "\t") %>%
  data.table::setnames(c("chr", "pos", "major", "minor")) %>%
  dplyr::select("chr", "pos") # don't need to keep which nucleotide ACGT is major/minor for hmm
map_pos = read.table(rdiff_file,
                     header = F, stringsAsFactors = F)
colnames(map_pos) = c("distM")
MEX_IDs = read.table(mex_file, stringsAsFactors = F, header = F)$V1
PARV_IDs = read.table(parv_file, stringsAsFactors = F, header = F)$V1
MAIZE_IDs = read.table(maize_file, stringsAsFactors = F, header = F)$V1

# helper function
# takes in individual counts file and samples 1 read per individual
# with either 0 0, 1 0 or 0 1 as values for n_major and n_minor at each SNP
sample1_read = function(ind_counts_file){
  counts = read.table(ind_counts_file, header = T, stringsAsFactors = F)
  tot = counts$n_major + counts$n_minor
  # if any reads are observed (tot > 0), sample 1 read. Otherwise count zero reads.
  new_n_minor = sapply(1:length(tot), function(x) ifelse(tot[x] == 0, 0, rbinom(n = 1,
                                            size = 1, # always sample 1 read
                                            prob = counts$n_minor[x]/tot[x])))
  new_n_major = ifelse(tot == 0, 0, 1 - new_n_minor)
  new_counts = data.frame(n_major = new_n_major, n_minor = new_n_minor)
  return(new_counts)
}

# start with 0 allele counts
counts = data.frame(n_major_mex = rep(0, nrow(SNPs)), 
                    n_minor_mex = rep(0, nrow(SNPs)),
                    n_major_parv = rep(0, nrow(SNPs)),
                    n_minor_parv = rep(0, nrow(SNPs)),
                    n_major_maize = rep(0, nrow(SNPs)),
                    n_minor_maize = rep(0, nrow(SNPs)))

# add mexicana allele counts
for (id in MEX_IDs){
  newCounts = sample1_read(
    ind_counts_file = paste0(DIR_COUNTS, "/", id, ".counts.txt"))
  colnames(newCounts) = c("n_major", "n_minor")
  # update count totals
  counts$n_minor_mex = counts$n_minor_mex + newCounts$n_minor
  counts$n_major_mex = counts$n_major_mex + newCounts$n_major
}

# add parviglumis allele counts
for (id in PARV_IDs){
  newCounts = sample1_read(
    ind_counts_file = paste0(DIR_COUNTS, "/", id, ".counts.txt"))
  colnames(newCounts) = c("n_major", "n_minor")
  # update count totals
  counts$n_minor_parv = counts$n_minor_parv + newCounts$n_minor
  counts$n_major_parv = counts$n_major_parv + newCounts$n_major
}

# add maize allele counts
for (id in MAIZE_IDs){
  newCounts = sample1_read(
    ind_counts_file = paste0(DIR_COUNTS, "/", id, ".counts.txt"))
  colnames(newCounts) = c("n_major", "n_minor")
  # update count totals
  counts$n_minor_maize = counts$n_minor_maize + newCounts$n_minor
  counts$n_major_maize = counts$n_major_maize + newCounts$n_major
}

# add SNP information to counts
d = bind_cols(SNPs, counts, map_pos)

# write output to be used by make_input_ancestry_hmm.R (with no headers)
options(scipen = 999)
write.table(d, output_file,
    row.names = F, col.names = T, quote = F, sep = " ") # prints column names even though final file for ancestry_hmm won't include them
warnings() # print any warnings

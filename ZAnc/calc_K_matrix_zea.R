#!/usr/bin/env Rscript
library(dplyr)

# this script takes in population ancestry frequencies for all sympatric pops 
# for maize or mexicana and saves the K matrix

# load variables from Snakefile
anc_rdata = snakemake@input[["pop_anc"]]
# anc_rdata = "local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/maize.pops.anc.RData"
k_functions = snakemake@input[["k_functions"]]
# k_functions = "ZAnc/k_matrix.R"
file_out = snakemake@output[["K"]]
# file_out = "ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/maize.K.RData"

# get functions
source(k_functions)

# load data
load(anc_rdata)

K = make_K_calcs(t(anc))

save(K, file = file_out)

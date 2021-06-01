#!/usr/bin/env Rscript
library(dplyr)

# this script takes in population ancestry frequencies for all sympatric pops
# for maize or mexicana and saves the K matrix

# load variables from Snakefile
anc_rdata = snakemake@input[["pop_anc"]]
# anc_rdata = "local_ancestry/results/ancestry_hmm/HILO_MAIZE55/K2/Ne10000_yesBoot/anc/maize.pops.anc.RData"
k_functions = snakemake@input[["k_functions"]]
# k_functions = "ZAnc/k_matrix.R"
file_out = snakemake@output[["K_matrix"]]
# file_out = "ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/maize.K.RData"
K_admix = snakemake@params[["K_admix"]]
# K_admix = 2

# get functions
source(k_functions)

# load data
load(anc_rdata)

# define mixing ancestries
ancestries = c("mexicana", "maize", "parv")[1:K_admix]

# calculate K matrix for each mixing ancestry
K = lapply(ancestries, function(a) make_K_calcs(t(anc[[a]])))

names(K) = ancestries

save(K, file = file_out)

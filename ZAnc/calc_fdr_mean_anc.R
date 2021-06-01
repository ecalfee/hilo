#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)

# this script calculates FDRs to show outlier loci

# load variables from Snakefile
# zea = "maize"
fdr_functions = snakemake@input[["fdr_functions"]]
# fdr_functions = "ZAnc/FDR.R"
sim_file = snakemake@input[["sim"]]
# sim_file = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, ".MVN.RData")
anc_file = snakemake@input[["anc"]]
# anc_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/K2/Ne10000_yesBoot/anc/", zea, ".pops.anc.RData")
fdr_out = snakemake@output[["fdr"]]
# fdr_out = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, ".meanAnc.fdr.RData")
K_admix = snakemake@params[["K_admix"]]
# K_admix = 2

# load data
source(fdr_functions)
load(anc_file)
load(sim_file)

# define mixing ancestries
ancestries = c("mexicana", "maize", "parv")[1:K_admix]

# calculate false discovery rate thresholds
FDRs = list()
for (a in ancestries){
  FDRs[[a]] = calc_FDR(d = anc_mean[[a]], 
                  s = MVN_mean[[a]], 
                  test_values = seq(0, 1, by = .0001))
  # what % of SNPs exceed these thresholds?
  FDRs[[a]]$n_SNPs = sapply(1:nrow(FDRs[[a]]), function(i) 
    ifelse(FDRs[[a]]$tail[i] == "high", 
           sum(anc_mean[[a]] > FDRs$threshold[i]),
           sum(anc_mean[[a]] < FDRs$threshold[i])))
  FDRs[[a]]$prop_SNPs = FDRs[[a]]$n_SNPs/length(anc_mean[[a]])
}

save(FDRs, file = fdr_out)
#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)

# this script calculates FDRs to show outlier loci

# load variables from Snakefile
# zea = "maize"
fdr_functions = snakemake@input[["fdr_functions"]]
# fdr_functions = "ZAnc/FDR.R"
sim_file = snakemake@input[["sim"]]
# sim_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".MVN.RData")
anc_file = snakemake@input[["anc"]]
# anc_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pops.anc.RData")
fdr_out = snakemake@output[["fdr"]]
# fdr_out = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".meanAnc.fdr.RData")

# load data
source(fdr_functions)
load(anc_file)
load(sim_file)

# calculate false discovery rate thresholds
FDRs = calc_FDR(d = anc_mean, 
                s = MVN_mean, 
                test_values = seq(0, 1, by = .0001))

# what % of SNPs exceed these thresholds?
FDRs$n_SNPs = sapply(1:nrow(FDRs), function(i) 
  ifelse(FDRs$tail[i] == "high", 
         sum(anc_mean > FDRs$thesholds[i]),
         sum(anc_mean < FDRs$thesholds[i])))
FDRs$prop_SNPs = FDRs$n_SNPs/length(anc_mean)

save(FDRs, file = fdr_out)
#!/usr/bin/env Rscript
library(dplyr)
# load variables from Snakefile
source(snakemake@input[["fdr_functions"]])
# source("ZAnc/FDR.R")
source(snakemake@input[["lm_functions"]])
# source("ZAnc/lm_env_function.R")
# zea = "maize"
meta_file = snakemake@input[["meta_pop"]]
# meta_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pop.meta.RData")
anc_file = snakemake@input[["anc"]]
# anc_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pops.anc.RData")
sim_file = snakemake@input[["sim"]]
# sim_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".MVN.RData")
fdr_out = snakemake@output[["fdr"]]
# fdr_out = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".lmElev.fdr.RData")
fit_out = snakemake@output[["fit"]]
# fit_out = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".lmElev.fit.RData")
sim_out = snakemake@output[["sim"]]
# sim_out = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".lmElev.sim.RData")

# load data
load(anc_file)
load(meta_file)
load(sim_file)


# fit original data
fits <- apply(anc, 
              1, function(x)
                simple_env_regression(ancFreq = x, 
                                      envWeights = meta_pops$ELEVATION/1000)) %>%
  t(.) %>%
  as.data.frame(.)

# fit simulated data
fits_sim <- apply(MVN_sim, 
              1, function(x)
                simple_env_regression(ancFreq = x, 
                                      envWeights = meta_pops$ELEVATION/1000)) %>%
  t(.) %>%
  as.data.frame(.)

# calculate false discovery rates
FDRs <- calc_FDR(d = fits$envWeights, s = fits_sim$envWeights, 
                 FDR_values = FDR_range, 
                             test_values = seq(min(fits$envWeights), 
                                               max(fits$envWeights), 
                                               by = .0001))

# what % of SNPs exceed these thresholds?
FDRs$n_SNPs = sapply(1:nrow(FDRs), function(i) 
  ifelse(FDRs$tail[i] == "high", 
         sum(fits$envWeights > FDRs$thesholds[i]),
         sum(fits$envWeights < FDRs$thesholds[i])))
FDRs$prop_SNPs = FDRs$n_SNPs/nrow(fits)

# save results
save(FDRs, file = fdr_out)
save(fits, file = fit_out)
save(fits_sim, file = sim_out)

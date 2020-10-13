#!/usr/bin/env Rscript

# this script calculates false discovery rate cutoffs for zAnc models
# based on enrichment in the data compared to the simulated MVN null model

library(dplyr)
# load variables from Snakefile
source(snakemake@input[["fdr_functions"]])
# source("ZAnc/FDR.R")
fdr_out = snakemake@output[["fdr"]]
# fdr_out = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".zAnc.fdr.RData")
sim_in = snakemake@output[["sim"]]
# sim_in = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".zAnc.sim.RData")
fit_in = snakemake@output[["fit"]]
# fit_in = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".zAnc.fit.RData")

# load data
load(sim_in)
load(fit_in)

# calculate false discovery rates
# use neutral MVN simulations, not theoretical chisq(df = 1) 
# to set significance at 5% FDR

# range of values to test go from top 10% of null simulations
# (all of test data would have to be above that threshold to meet 10% FDR)
# to maximum value in test data
print("ranges tested for FDR:")
print("zElev")
print(quantile(fits_sim$diff_elev_ll2, 1 - max(FDR_range)))
print(max(fits$diff_elev_ll2))
print("zAll")
print(quantile(fits_sim$diff_all_ll2, 1 - max(FDR_range)))
print(max(fits$diff_all_ll2))
print("zTz")
print(quantile(fits_sim$zTz, 1 - max(FDR_range)))
print(max(fits$zTz))

# calculate false-discovery rate -- unusually high values indicate statistical significance
FDRs <- bind_rows(calc_FDR_high(d = fits$diff_elev_ll2, s = fits_sim$diff_elev_ll2, 
                                FDR_values = FDR_range, 
                                test_values = seq(quantile(fits_sim$diff_elev_ll2, 1 - max(FDR_range)), 
                                                  max(fits$diff_elev_ll2), 
                                                  by = .0001)) %>%
                    dplyr::mutate(model = "elev"),
                  calc_FDR_high(d = fits$diff_all_ll2, s = fits_sim$diff_all_ll2, 
                                FDR_values = FDR_range, 
                                test_values = seq(quantile(fits_sim$diff_all_ll2, 1 - max(FDR_range)), 
                                                  max(fits$diff_all_ll2), 
                                                  by = .0001)) %>%
                    dplyr::mutate(model = "all"),
                  calc_FDR_high(d = fits$zTz, s = fits_sim$zTz, 
                                FDR_values = FDR_range, 
                                test_values = seq(quantile(fits_sim$zTz, 1 - max(FDR_range)), 
                                                  max(fits$zTz), 
                                                  by = .0001)) %>%
                    dplyr::mutate(model = "zTz"))

# what % of SNPs exceed these thresholds?
FDRs$n_SNPs = sapply(1:nrow(FDRs), function(i) 
  sum(fits[ , paste0("diff_", FDRs$model[i], "_ll2")] > FDRs$thesholds[i]))
FDRs$prop_SNPs = FDRs$n_SNPs/nrow(fits)

# save results
save(FDRs, file = fdr_out)
#!/usr/bin/env Rscript

# this script fits zAnc models to data and simulated null:
# zElev (elevation), zAll (all pops selected + or -), zTz (general ancestry outlier test)

library(dplyr)
# load variables from Snakefile
source(snakemake@input[["zAnc_functions"]])
# source("ZAnc/other_functions.R")
# zea = "maize"
meta_file = snakemake@input[["meta_pop"]]
# meta_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pop.meta.RData")
anc_file = snakemake@input[["anc"]]
# anc_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pops.anc.RData")
sim_file = snakemake@input[["sim"]]
# sim_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".MVN.RData")
k_file = snakemake@input[["K"]]
# k_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".K.RData")
sim_out = snakemake@output[["sim"]]
# sim_out = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".zAnc.sim.RData")
fit_out = snakemake@output[["fit"]]
# fit_out = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".zAnc.fit.RData")

# load data
load(anc_file)
load(meta_file)
load(sim_file)
load(k_file) # k matrix

# calculate elevation in km, mean-centered across pops
meta_pops$elev_km <- meta_pops$ELEVATION/1000
meta_pops$c_elev_km <- meta_pops$elev_km - mean(meta_pops$elev_km)

# solve for inverse of K matrix
invK = solve(K$K)
# solve for determinant of K matrix
detK = det(K$K)

# function to fit all selection scenarios
fit_zAnc <- function(anc, alpha, invK, detK, c_elev_km){
  fits <- list(null = fit_null(anc = anc, 
                                 alpha = alpha, invK = invK, detK = detK) %>%
                 data.table::setnames(., paste0(colnames(.), ".null")), 
                 allpop = fit_sel(anc = anc, 
                                  alpha = alpha, invK = invK, detK = detK, 
                                  X = as.matrix(rep(1, length(alpha))), 
                                  b_names = "b") %>%
                 data.table::setnames(., paste0(colnames(.), ".all")), 
                 elevation = fit_sel(anc = anc, 
                                     alpha = alpha, invK = invK, detK = detK, 
                                     X = as.matrix(cbind(intercept = 1, c_elev_km)), 
                                     b_names = c("b0", "b1")) %>%
                 data.table::setnames(., paste0(colnames(.), ".elev"))) %>%
    do.call(bind_cols, .) %>% 
    # using likelihood ratio, compare:
    # slope + intercept model (elevation) vs. just intercept model (allpop) and, separately,
    # just intercept model (allpop) vs. no intercept (null)
    dplyr::mutate(diff_elev_ll2 = -2*(ll.all - ll.elev), # calculate 2*ln(likelihood ratio)
                  diff_all_ll2 = -2*(ll.null - ll.all),
                  zTz = RSS.null) # zTz statistic is just the residual sum of squares under the null
  return(fits) 
}

# fit original data
fits <- fit_zAnc(anc = data.frame(anc), alpha = K$alpha, invK = invK, detK = detK, c_elev_km = meta_pops$c_elev_km)

# fit simulated data
fits_sim <- fit_zAnc(anc = data.frame(MVN_sim), alpha = K$alpha, invK = invK, detK = detK, c_elev_km = meta_pops$c_elev_km)

# save results
save(fits, file = fit_out)
save(fits_sim, file = sim_out)
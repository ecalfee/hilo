#!/usr/bin/env Rscript
library(dplyr)
library(MASS)

# this script simulates mexicana ancestry frequencies
# based on a ~MVN(alpha, K) with truncation to range [0,1]

# load variables from Snakefile
# zea = "maize"
K_file = snakemake@input[["K"]]
# K_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".K.RData")
meta_file = snakemake@input[["meta_pop"]]
# meta_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pop.meta.RData")
file_out = snakemake@output[["mvn"]]
# file_out = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".MVN.RData")
n_sim = as.integer(snakemake@params[["n_sim"]])
# n_sim = 10000

# load data
load(K_file)
load(meta_file)

# create some MVN data
set.seed(100)
MVN_sim_untruncated = mvrnorm(n = n_sim, 
                              mu = K$alpha,
                              Sigma = K$K)

# proportion truncated
print("Truncated too low:")
apply(MVN_sim_untruncated, 2, function(x) sum(x < 0))/n_sim
print("Truncated too high:")
apply(MVN_sim_untruncated, 2, function(x) sum(x > 1))/n_sim

# truncate simulated allele freqs within [0,1] range
truncate01 <- function(x){
  t <- x
  t[t < 0] <- 0
  t[t > 1] <- 1
  return(t)
}
MVN_sim = truncate01(MVN_sim_untruncated)

# calculate mean across individuals
MVN_mean = (MVN_sim %*% meta_pops$n_local_ancestry)/sum(meta_pops$n_local_ancestry)

save(MVN_sim, MVN_mean, file = file_out)

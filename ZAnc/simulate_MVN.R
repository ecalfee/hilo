#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(MASS)

# this script simulates mexicana ancestry frequencies
# based on a ~MVN(alpha, K) with truncation to range [0,1]

# load variables from Snakefile
# zea = "maize"
K_file = snakemake@input[["K"]]
# K_file = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, ".K.RData")
meta_file = snakemake@input[["meta_pop"]]
# meta_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/K2/Ne10000_yesBoot/anc/", zea, ".pop.meta.RData")
file_out = snakemake@output[["sim"]]
# file_out = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, ".MVN.RData")
n_sim = as.integer(snakemake@params[["n_sim"]])
# n_sim = 10000
K_admix = snakemake@params[["K_admix"]]
# K_admix = 2
truncated_txt = snakemake@output[["truncated_txt"]]
# truncated_txt = paste0("ZAnc/results/HILO_MAIZE55/K2/Ne10000_yesBoot/", zea, ".MVN.truncated.stats.txt")

# load data
load(K_file)
load(meta_file)

# define mixing ancestries
ancestries = c("mexicana", "maize", "parv")[1:K_admix]

# create some MVN data
set.seed(100)
MVN_sim_untruncated = list()
for (a in ancestries){
  MVN_sim_untruncated[[a]] = mvrnorm(n = n_sim,
                                     mu = K[[a]]$alpha,
                                     Sigma = K[[a]]$K)
}

# what proportion of simulated ancestry frequencies per population fall outside of the range [0,1]?
trunc_low = do.call(bind_rows, lapply(ancestries, function(a)
  apply(MVN_sim_untruncated[[a]], 2, function(x) sum(x < 0))/n_sim)) %>%
  dplyr::mutate(ancestry = ancestries) %>%
  tidyr::pivot_longer(data = ., 
                      cols = starts_with("pop"),
                      names_to = "pop",
                      values_to = "too_low")
trunc_high = do.call(bind_rows, lapply(ancestries, function(a)
  apply(MVN_sim_untruncated[[a]], 2, function(x) sum(x > 1))/n_sim)) %>%
  dplyr::mutate(ancestry = ancestries) %>%
  tidyr::pivot_longer(data = ., 
                      cols = starts_with("pop"),
                      names_to = "pop",
                      values_to = "too_high")
bind_rows(trunc_low, trunc_high) %>%
  write.table(x = ., file = truncated_txt, 
              col.names = T, row.names = F,
              sep = "\t", quote = F)

# truncate simulated allele freqs to within [0,1] range
truncate01 <- function(x){
  t <- x
  t[t < 0] <- 0
  t[t > 1] <- 1
  return(t)
}

MVN_sim = list()
MVN_mean = list()

for (a in ancestries){
  # truncated simulation values
  MVN_sim[[a]] = truncate01(MVN_sim_untruncated[[a]])
  
  # calculate mean across individuals
  MVN_mean[[a]] = (MVN_sim[[a]] %*% meta_pops$n_local_ancestry)/sum(meta_pops$n_local_ancestry)
}

# save results
save(MVN_sim, MVN_mean, file = file_out)

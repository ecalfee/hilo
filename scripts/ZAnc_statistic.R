library(dplyr)
library(ggplot2)
# helper file with some useful functions
source("../../covAncestry/forqs_sim/k_matrix.R") # import useful functions

dir_in = "../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/"
dir_post = paste0(dir_in, "output_noBoot/")
# admixture proportions per pop
alphas <- read.table("../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/input/globalAdmixtureByPopN.txt",
                     stringsAsFactors = F, header = F)
colnames(alphas) <- c("popN", "alpha_maize", "alpha_mex")
# metadata for each individual, including pop association
meta <- read.table("../data/pass1_ids.txt", stringsAsFactors = F, 
                   header = T, sep = "\t") %>%
  filter(est_coverage >= 0.05) %>%
  left_join(., alphas, by = c("popN"))
# individual by individual K matrix
included_inds = meta %>%
  filter(alpha_maize > 0 & alpha_mex >0)
sample_pos = function(nth, ID, path) {# take every nth SNP & summarize across posterior of all 3 genotypes 
  # to get estimate of mex. ancestry at each of those SNPs for a HILO individual ID
  as.matrix(read.table(paste0(path, "HILO", ID, ".posterior"), 
                       header = T, 
                       stringsAsFactors = F)[c(T, rep(F, nth-1)), 3:5]) %*% 
    c(0, .5, 1) 
}
# get (mexicana) ancestry in the form 
# 1 sample per row and 1 marker per column
#listOfIds = included_inds$n
listOfIds = included_inds$n[included_inds$popN == 366]
anc = t(do.call(cbind, 
          lapply(listOfIds, 
                 function(i) sample_pos(nth = 50, ID = i, path = dir_post))))
alpha_ind = apply(anc, 1, mean) # mean mexicana ancestry across all markers per individual
# make separate script that creates one ancestry file per population (all SNPs)

# maybe make list of pops, where each element of the list is a vector of pop IDs
# then I can convert to a list of 'anc' matrices per pop, which can be rbind to get ind ancestry
# or first take mean within pop then combine across pops into 1 matrix
# short goal 1: plot of individual K matrix for all maize, ordered by location
# all mex, ordered by location;
# and then combined maize/mex samples, ordered by location
# later goal: low vs. high recombination regions, same K matrices
# short goal 2: calculate population-level 
#anc_pop = tapply(anc$, )
mex = as.matrix(post_mex[, 3:5]) %*% c(0, .5, 1) 
# population by population K matrix


# create K matrices and calc alpha and inverse cholesky
K = lapply(sims, function(i) 
  make_K_matrix(get_freqs(
    paste0(pathAnc, prefix, ".sim", i, "/", nChr, "/population_", gen, "_", pops, ".frq")),
    ploidy = nHap)) # the population ploidy is the # of haplotypes with ancestry counts per population

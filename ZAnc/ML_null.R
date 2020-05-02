#!/usr/bin/env Rscript
# this script calculates the likelihood under the null model
# ancestry ~ MVN(mu, K), where mu = alpha;
# no parameters to estimate
# returns log likelihood and AICc

# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# functions
source("k_matrix.R")

# load data 
PREFIX="pass2_alloMAIZE"

# metadata for each individual, including pop association
pop_elev <- read.table("../data/riplasm/gps_and_elevation_for_sample_sites.txt",
                       stringsAsFactors = F, header = T, sep = "\t") %>%
  dplyr::select(., popN, ELEVATION, LAT, LON)
hilo <- read.table("../samples/hilo_meta.txt", stringsAsFactors = F, header = T, sep = "\t")
# admixture proportions per pop
meta.pops <- read.table(paste0("../global_ancestry/results/NGSAdmix/", PREFIX, "/globalAdmixtureIncludedByPopN.txt"),
                        stringsAsFactors = F, header = T) %>%
  left_join(., pop_elev, by = c("popN")) %>%
  left_join(., unique(dplyr::select(hilo, c("popN", "zea", "symp_allo", "RI_ACCESSION", "GEOCTY", "LOCALITY"))), by = c("popN")) %>%
  mutate(pop = paste0("pop", popN))

# ancestry calls
LOCAL_ANC_SUBDIR="output_noBoot"
dir_in = paste0("../local_ancestry/results/ancestry_hmm/", PREFIX)
dir_anc = file.path(dir_in, LOCAL_ANC_SUBDIR, "anc")
# read in population ancestry input files from calc_genomewide_pop_anc_freq.R
pop_anc_list = lapply(meta.pops$popN, function(pop) read.table(paste0(dir_anc, "/pop", pop, ".anc.freq"), 
                                                               stringsAsFactors = F))

# combine population ancestry frequencies into a matrix
# where rows are populations and columns are snps
all_anc = do.call(cbind, pop_anc_list)
colnames(all_anc) <- meta.pops$pop
rm(pop_anc_list)

# get mean ancestry across individuals sampled (NOT mean of pops)
maize_pops <- unique(meta.pops$pop[meta.pops$zea == "maize"])
#mexicana_pops <- unique(meta.pops$pop[meta.pops$zea == "mexicana"])
maize_anc <- all_anc[ , maize_pops]
#mexicana_anc <- all_anc[ , mexicana_pops]

# calculations
# ancestry mean and covariance matrix
#zAnc_mexicana = make_K_calcs(t(mexicana_anc))
zAnc_maize = make_K_calcs(t(maize_anc))
alpha = zAnc_maize$alpha
K = zAnc_maize$K
invK = solve(K)
detK = det(K) # determinant of K matrix

# estimate log likelihood under this MVN model
logliks <- sapply(1:nrow(maize_anc), function(i)
  ll_mvn(t(maize_anc[i, ]), # make into a column vector
         mu = alpha,
         detK = detK, 
         invK = invK))

rss <- sapply(1:nrow(maize_anc), function(i)
  RSS(ancFreq = t(maize_anc[i, ]),
         mu = alpha,
         invK = invK))

fits <- data.frame(
  RSS = rss,
  ll = logliks,
  k = 0,
  n = length(alpha)) %>%
  mutate(AICc = AICc(ll = ll, k = k, n = n))

write.table(fits,
            paste0("results/models/", PREFIX, "/maize/ML_null.txt"),
            sep = "\t",
            quote = F,
            col.names = T, 
            row.names = F)

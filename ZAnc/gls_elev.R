#!/usr/bin/env Rscript
# this script calculates generalized least squares results at all loci across
# maize genome for a linear model with ancestry ~ b0 + bElev*elev_c
# where elev_c is the population elevation (in km), centered at the mean elev across pops
# returns b0 and bElev estimates for each snp

# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

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

# genomic sites information
# SNPs with ancestry calls
dir_sites = paste0("../local_ancestry/results/thinnedSNPs/", PREFIX)
sites <- do.call(rbind, 
                 lapply(1:10, function(i)
                   read.table(paste0(dir_sites, "/chr", i, ".var.sites"),
                              header = F, stringsAsFactors = F)))
colnames(sites) <- c("chr", "pos", "major", "minor")
# get inversions:
inv = read.table("../data/refMaize/inversions/knownInv_v4_coord.txt",
                 stringsAsFactors = F, header = T)
colnames(inv) = c("ID", "chr", "start", "end", "length")
excl_inv_buffer = 1000000 # exclude 1 megabase around inversion
inv$excl_start = inv$start - excl_inv_buffer
inv$excl_end = inv$end + excl_inv_buffer
# which sites are within the inversion?
sites <- mutate(sites, inv4m_1Mb_buffer = (chr == inv$chr[inv$ID=="inv4m"] & 
                                             pos >= inv$excl_start[inv$ID=="inv4m"] & 
                                             pos <= inv$excl_end[inv$ID == "inv4m"]),
                inv4m = (chr == inv$chr[inv$ID=="inv4m"] & 
                           pos >= inv$start[inv$ID=="inv4m"] & 
                           pos <= inv$end[inv$ID == "inv4m"])) 


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
#zAnc_mexicana = make_K_calcs(t(mexicana_anc))
zAnc_maize = make_K_calcs(t(maize_anc))
alpha = zAnc_maize$alpha
K = zAnc_maize$K
elev = sapply(names(alpha), function(p)
  meta.pops$ELEVATION[meta.pops$pop == p]/1000)
elev_c = elev - mean(elev)

# what I really want to do is not a logistic, but add slope, then truncation at 0, 1
# but logistic might be a good alternative. apply logistic before or after noise?
X = cbind(intercept = 1, elev_c) %>%
  as.matrix(.)
invK = solve(K)

# generalized least squares model for population ancestries centered at pop mean ancestry:
# (y - alpha)
# take out part1 of the generalized least squares model to not repeat calcs:
part1 = solve(t(X) %*% invK %*% X) %*% t(X) %*% invK
# now apply to each observed vector of population ancestries
betas = t(apply(maize_anc, 1, function(y) 
  part1 %*% (y - alpha))) %>%
  as.data.frame(.) %>%
  data.table::setnames(c("b0", "bElev"))

write.table(betas,
            paste0("results/models/", PREFIX, "/maize/gls_elevation.txt"),
            sep = "\t",
            quote = F,
            col.names = T, 
            row.names = F)

#zAnc3_mex <- t(apply(mexicana_anc, 1, function(l) zAnc(ancFreq = l, 
#                                                       invL = zAnc_mexicana$InvL, 
#                                                       alpha = zAnc_mexicana$alpha)))
#colnames(zAnc3_mex) <- colnames(mexicana_anc)

# simple simulations of environmental associations + MVN noise
library(dplyr)
library(tidyr)
library(ggplot2)

# read in data (or could pick arbitrary alpha's and K matrices)
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

maize_pops <- unique(meta.pops$pop[meta.pops$zea == "maize"])


K <- read.table("results/zAnc_K_matrix_maize.txt", sep = "\t", header = T)
alpha <- t(read.table("results/zAnc_alpha_maize.txt", sep = "\t", header = T))

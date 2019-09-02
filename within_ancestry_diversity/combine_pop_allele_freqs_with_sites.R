#!/usr/bin/env Rscript
# this script loads population allele frequencies
# for a pop, and combines them with the sites file,
# producing a new file with a frequency for every site (or NA if no data)
# and number of individuals with data for each site

library(dplyr)

# to run:
# Rscript combine_pop_allele_freqs_with_sites.R parv results/allele_freq/pass2_alloMAIZE/parv/region_0 ../variant_sites/results/pass2_alloMAIZE/region_0.var.sites

# arguments
args = commandArgs(trailingOnly=TRUE)
# population number
POP = args[1]
# paths to input and output directories
path_allele_freqs = args[2]
# sites file
sites_file = args[3]

sites0 <- read.table(sites_file,
                    stringsAsFactors = F, header = F) %>%
  data.table::setnames(c("scaffold", "pos", "major", "minor"))

allele_freq <- read.table(paste0(path_allele_freqs, ".mafs.gz"),
                          stringsAsFactors = F, header = T) %>%
  left_join(sites0, ., by = c("scaffold"="chromo", "pos"="position",
                              "major", "minor"))

# separate out allele freqs and write to file
f <- allele_freq %>%
  mutate(p = ifelse(phat > 1, 1, phat)) %>% # gets rid of weird angsd rounding error
  dplyr::select(p) %>%
  data.table::setnames(POP)
write.table(f, paste0(path_allele_freqs, ".freqs.txt"),
            col.names = T, row.names = F, quote = F)

# separate out number of individuals with data and write to file
n <- allele_freq %>%
  mutate(n = ifelse(is.na(nInd), 0, nInd)) %>% # NA means no individuals with data
  dplyr::select(n) %>%
  data.table::setnames(POP)
write.table(n, paste0(path_allele_freqs, ".nInd"),
            col.names = T, row.names = F, quote = F)

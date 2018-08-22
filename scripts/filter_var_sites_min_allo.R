library(dplyr)
# this script filters var sites before an LD filter
# for ultimately inputting into ancestry_hmm for local ancestry calling

minInd = 4 # minimum number of individuals per population to include a site
D = .3 # minimum difference in estimated allele frequency between 
# allopatric maize and allopatric mexicana (goal = to enrich for informative sites)

# prefix/main directory for data
prefix = "../data/geno_lik/merged_pass1_all_alloMaize4Low_16/allVar/"
region = "region_9"
# get allopatric MAF data for some region
maize = read.table(gzfile(paste0(prefix, "maize.allo.4Low16", "/", 
                                 region, ".mafs.gz")), header=T, stringsAsFactors = F) %>%
  rename(.,  MAF_maize = knownEM) %>%
  rename(., nInd_maize = nInd) %>%
  select(., -ref) # remove reference

mex = read.table(gzfile(paste0(prefix, "mexicana.allo.withXochi35", "/",
                               region, ".mafs.gz")), header=T, stringsAsFactors = F) %>%
  rename(.,  MAF_mex = knownEM) %>%
  rename(., nInd_mex = nInd) %>%
  select(., -ref) # remove reference

#d2 = read.table(gzfile(paste0(prefix, "whole_genome_pruned_every_1000.beagle.gz")), header=T, stringsAsFactors = F)
#sum(duplicated(d2[ , 1:4])) # no duplicated sites by chance in the every 1000 SNPs file

d = read.table(paste0(prefix, region, ".var.sites"), header=F, stringsAsFactors = F) %>%
  setNames(., colnames(maize)[1:4]) %>%
  left_join(., maize, by = colnames(maize)[1:4]) %>%
  left_join(., mex, by = colnames(maize)[1:4]) %>%
  filter(., !duplicated(.)) # temporary problem with duplicated positions -- I'll fix upstream later

# try different filtering and look at results
tot = sum(complete.cases(d)) # number of sites with data from maize and mex (largest possible total)
d %>% # what % of sites will be left after some minor filtering for more informative sites?
  filter(., nInd_maize >= minInd & nInd_mex >= minInd) %>%
  filter(., abs(MAF_maize - MAF_mex) >= D) %>%
  nrow(.)/tot
# at a filter of at least 4 individuals per allopatric population,
# around .2 allele freq. difference or .3 is ~30% and ~17% of sites remaining post-filtering
# D = .4 drops the number of sites down to 1% 
# so my estimate is D = .3 will be best for ancestry calling to still get high enough
# density of SNPs but enriching for more informative SNPs .. I will still need to do some 
# filtering for SNPs in very high LD (too close cM distance)


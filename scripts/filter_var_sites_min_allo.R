#!/usr/bin/env Rscript

#R script 2 filter variant sites for minimum allopatric # individuals,
#MAF difference, and an LD threshold
print(getwd()) # print current directory
args = commandArgs(trailingOnly=TRUE)
# args (in order): region, minInd, D, allo_maize, allo_mex, dir_gl, dir_sites, dir_out
# to run:
# Rscript filter_var_sites_min_allo.R 9 4 0.3 maize.allo.4Low16 mexicana.allo.withXochi35 (args cont. next lines..)
# ../data/geno_lik/merged_pass1_all_alloMaize4Low_16/allVar/
# ../data/var_sites/merged_pass1_all_alloMaize4Low_16/
# ../data/var_sites/merged_pass1_all_alloMaize4Low_16/filteredAlloMAFnInd

# If I wanted to convert this to a python script I could use the pandas package, e.g.
# pd.merge(df1, df2, on=["chromo", "position", "major", "minor"],
# suffices=["_mex", "_maize"], how="inner") # or outer join for other applications
library(dplyr)
# this script filters var sites before an LD filter
# for ultimately inputting into ancestry_hmm for local ancestry calling
#region = 9 # region of the genome
region = as.integer(args[1])
#minInd = 4 # minimum number of individuals per population to include a site
minInd = as.integer(args[2])
#D = 0.3
D = as.numeric(args[3])
# minimum difference in estimated allele frequency between
# allopatric maize and allopatric mexicana (goal = to enrich for informative sites)

# prefix/main directory for data
#allo_maize = "maize.allo.4Low16"
allo_maize = args[4]
#allo_mex = "mexicana.allo.withXochi35"
allo_mex = args[5]
#dir_maf = "../data/var_sites/merged_pass1_all_alloMaize4Low_16/"
dir_maf = args[6]
#dir_sites = "../data/var_sites/merged_pass1_all_alloMaize4Low_16/"
dir_sites = args[7]
#dir_out = "../data/var_sites/merged_pass1_all_alloMaize4Low_16/filteredAlloMAFnInd"
dir_out = args[8]

# function to get variant sites and allopatric maize and mexicana frequencies in one dataframe
getSNPs = function(region, dir_gl, dir_sites, allo_maize, allo_mex, D, minInd){
  # get allopatric MAF data for some region
  maize = read.table(gzfile(paste0(dir_maf, allo_maize, "/",
                                   "region_", region, ".mafs.gz")), header=T, stringsAsFactors = F) %>%
    rename(.,  MAF_maize = knownEM) %>%
    rename(., nInd_maize = nInd) %>%
    select(., -ref) # remove reference

  mex = read.table(gzfile(paste0(dir_maf, allo_mex, "/",
                                 "region_", region, ".mafs.gz")), header=T, stringsAsFactors = F) %>%
    rename(.,  MAF_mex = knownEM) %>%
    rename(., nInd_mex = nInd) %>%
    select(., -ref) # remove reference

  #d2 = read.table(gzfile("../data/geno_lik/merged_pass1_all_alloMaize4Low_16/allVar/whole_genome_pruned_every_1000.beagle.gz"), header=T, stringsAsFactors = F)
  #sum(duplicated(d2[ , 1:4])) # no duplicated sites by chance in the every 1000 SNPs file

  d = read.table(paste0(dir_sites, "/", "region_", region, ".var.sites"),
                 header=F, stringsAsFactors = F) %>%
    setNames(., colnames(maize)[1:4]) %>%
    left_join(., maize, by = colnames(maize)[1:4]) %>%
    left_join(., mex, by = colnames(maize)[1:4]) %>%
    filter(., !duplicated(.)) %>% # temporary problem with duplicated positions -- I'll fix upstream later
    filter(., nInd_maize >= minInd & nInd_mex >= minInd) %>% # thin SNPs based on min # ind's in allopatric pops
    filter(., abs(MAF_maize - MAF_mex) >= D) # and minimum allele freq. difference maize-mex
  return(d)
}


# try different filtering and look at results
# what % of sites will be left after some minor filtering for more informative sites?
#tot = sum(complete.cases(d)) # number of sites with data from maize and mex (largest possible total)
#nrow(d_thin)/tot
# at a filter of at least 4 individuals per allopatric population,
# around .2 allele freq. difference or .3 is ~30% and ~17% of sites remaining post-filtering
# D = .4 drops the number of sites down to 1%
# so my estimate is D = .3 will be best for ancestry calling to still get high enough
# density of SNPs but enriching for more informative SNPs .. I will still need to do some
# filtering for SNPs in very high LD (too close cM distance)
#sum(keep)
#hist(d_thin$nInd_maize)
#hist(d_thin$nInd_mex)
#table(d_thin$nInd_mex)

# filter SNPs for informative sites
d = getSNPs(region = region, dir_maf = dir_maf,
            dir_sites = dir_sites, allo_maize = allo_maize,
            allo_mex = allo_mex, D = D, minInd = minInd)
# write filtered SNPs to a new sites file
write.table(file = paste0(dir_out, "/region_", region, ".var.sites"),
            x = d[ , c("chromo", "position", "major", "minor")],
            quote = F, row.names = F, col.names = F, sep = "\t")

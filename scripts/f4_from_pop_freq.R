#!/usr/bin/env Rscript
library(dplyr)
library(bootstrap)
# this script runs an f4 test of admixture
# Input:
# It takes in the names of 4 populations
# to run F4(pop1, pop2, pop3, pop4),
# a directory for where to find their angsd-output
# allele freq. files for a set of regions, 
# e.g. dir_in/pop_name/region_i.mafs.gz
# Output: 
# A results file with F4 test and leave-one-out jackknife
# bias and s.e.
# As well as D (normalized F4) calculated as in doAbbaBaba2 in Angsd
# And both F3 statistics, testing for admixture in pop2 and pop3

# A regions file with 3 columns: 
# (1) region number
# (3) # variable sites across 4 pops in region
# (2) sum of (x1-x2)(x3-x4)
# (4) sum of (x1+x2-2*x1*x2)(x3+x4-2*x3*x4)
# (5) sum of (x2 - x1)*(x2 - x4) #F3(pop2; pop1, pop4)
# (6) sum of (x3 - x1)*(x3 - x4) #F3(pop3; pop1, pop4)

# to run:
# Rscript f4_from_pop_freq.R maize.allo.4Low16 maize.symp mexicana.symp mexicana.allo ../data/geno_lik/merged_pass1_all_alloMaize4Low_16/allVar_depthFilt/popFreqs

## NOTE (!) not sure why but I'm getting many values 
# 1.001002 > freq greater than 1 (rounding error in angsd?)

# arguments
regions <- 0:425 # regions to analyze (5Mb each)
#regions <- 20:21

args = commandArgs(trailingOnly=TRUE)
#args = c("maize.allo.4Low16", "maize.symp", 
#         "mexicana.symp", "mexicana.allo", 
#         "../data/geno_lik/merged_pass1_all_alloMaize4Low_16/allVar_depthFilt/popFreqs")
# population number
pops = args[1:4]
dir = args[5]
dir_output = paste0(dir, "/f4")
if (file.exists(dir) & !file.exists(dir_output)){ # make sure output directory exists first (or create)!
  dir.create(file.path(dir_output), recursive = T)
}

file_out = paste0(dir_output, "/", pops[1], "_", 
                  pops[2], "_", pops[3], "_", pops[4])

# helper function to calculate f4 for genomic region i
calc_f4_by_region <- function(i){ 
  # read in population allele freq files
  pop1 <- read.table(paste0(dir, "/", pops[1], "/region_",
                            i, ".mafs.gz"), header = T)[ , c(1:4,6)] %>%
    rename(x1 = phat)
  pop2 <- read.table(paste0(dir, "/", pops[2], "/region_",
                            i, ".mafs.gz"), header = T)[ , c(1:4,6)] %>%
    rename(x2 = phat)
  pop3 <- read.table(paste0(dir, "/", pops[3], "/region_",
                            i, ".mafs.gz"), header = T)[ , c(1:4,6)] %>%
    rename(x3 = phat)
  pop4 <- read.table(paste0(dir, "/", pops[4], "/region_",
                            i, ".mafs.gz"), header = T)[ , c(1:4,6)] %>%
    rename(x4 = phat)
  
  # join pop freqs into one dataframe
  freqs <- pop1 %>%
    inner_join(., pop2, 
               by = c("chromo", "position", "major", "minor")) %>%
    inner_join(., pop3,
               by = c("chromo", "position", "major", "minor")) %>%
    inner_join(., pop4,
               by = c("chromo", "position", "major", "minor")) %>%
    filter(., x1 + x2 + x3 + x4 > 0 &
             x1 + x2 + x3 + x4 < 4) %>% # filter out invariant sites
    mutate(., f4 = (x1 - x2)*(x3 - x4)) %>%
    mutate(., D_denom = (x1 + x2 - 2*x1*x2)*(x3 + x4 - 2*x3*x4)) # denominator to normalize D statistic
  
  #F3(pop2; pop1, pop4) -- testing for admixture in pop2
  freqs_f3_2 <- freqs %>%
    filter(., x1 + x2 + x4 > 0 & x1 + x2 + x4 < 3) %>%
    mutate(., f3 = (x2 - x1)*(x2 - x4))
  
  #F3(pop3; pop1, pop4) -- testing for admixture in pop3
  freqs_f3_3 <- freqs %>%
    filter(., x1 + x3 + x4 > 0 & x1 + x3 + x4 < 3) %>%
    mutate(., f3 = (x3 - x1)*(x3 - x4))
  
  # calculate statistics
  n <- nrow(freqs) # number of loci included
  f4_sum <- sum(freqs$f4)
  D_denom_sum <- sum(freqs$D_denom)
  f3_sum2 <- sum(freqs_f3_2$f3)
  f3_sum3 <- sum(freqs_f3_3$f3)
  
  return(data.frame(region = i, n, f4_sum, D_denom_sum,
                    f3_sum2, f3_sum3))
}

# calculate and combine results for all reagions
d <- do.call(rbind,
             lapply(regions, function(i)
               calc_f4_by_region(i = i)))

# write output file for regions:
write.table(d, 
            paste0(file_out, ".regions"), 
            col.names = T, row.names = F, quote = F, sep = "\t")

# compute jackknife across regions
calc_f4 <- function(x, df){ # x is a vector of included rows
  sum(df[x, "f4_sum"])/sum(df[x, "n"])
}
jack_f4 <- jackknife(1:nrow(d), calc_f4, d)
f4 <- calc_f4(1:nrow(d), d) # with all data

# function to calculate normalized f4
# equivalent to D from doAbbaBaba2 in angsd
calc_D <- function(x, df){
  sum(df[x, "f4_sum"])/sum(df[x, "D_denom_sum"])
}
jack_D <- jackknife(1:nrow(d), calc_D, d)
D <- calc_D(1:nrow(d), d)

# both possible f3 tests:
calc_f3 <- function(x, df, f3_pop){ # x is a vector of included rows
  sum(df[x, f3_pop])/sum(df[x, "n"])
}
jack_f3 <- lapply(2:3, function(i)
                  jackknife(1:nrow(d), calc_f3, d, paste0("f3_sum", i)))
f3 <- lapply(2:3, function(i)
  calc_f3(1:nrow(d), d, paste0("f3_sum", i))) # with all data

# write output file for summary of jackknife
write.table(data.frame(pop1 = pops[1],
                       pop2 = pops[2],
                       pop3 = pops[3],
                       pop4 = pops[4],
                       f4 = f4,
                       f4.jack.se = jack_f4$jack.se,
                       f4.jack.bias = jack_f4$jack.bias,
                       D = D,
                       D.jack.se = jack_D$jack.se,
                       D.jack.bias = jack_D$jack.bias,
                       f3.pop2 = f3[[1]],
                       f3.pop2.jack.se = jack_f3[[1]]$jack.se,
                       f3.pop2.jack.bias = jack_f3[[1]]$jack.bias,
                       f3.pop3 = f3[[2]],
                       f3.pop3.jack.se = jack_f3[[2]]$jack.se,
                       f3.pop3.jack.bias = jack_f3[[2]]$jack.bias),
            paste0(file_out, ".summary"), 
            col.names = T, row.names = F, quote = F, sep = "\t")


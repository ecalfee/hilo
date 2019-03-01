#!/usr/bin/env Rscript

# this script takes as input a set of prefix name
# and looks for results/prefix/pop.anc.list file
# with paths to all the population ancestry files to include
# makes output file with ZAnc score for every SNP included
# results/prefix/pop.ZAnc

# to run: Rscript ./calc_ZAnc_from_pop_freq.R pass1_maize

# arguments
args = commandArgs(trailingOnly=TRUE)
# prefix
prefix = args[1]

# helper file with some useful functions
source("../../covAncestry/forqs_sim/k_matrix.R") # import useful functions

# helper function for calculating statistic:
make_calcs = function(pop_anc){
  # calculate mean population frequency across all snps
  pop_alpha = apply(pop_anc, 1, mean)
  # calculate K matrix for populations
  pop_K = calcK(ancFreqMatrix = pop_anc, alpha = pop_alpha)
  pop_InvL = calcInvL(pop_K)
  mean_anc = apply(pop_anc, 2, mean) # mean across populations
  # for each locus, calculate ZAnc statistic
  pop_ZAnc = apply(pop_anc, 2, function(l) ZAnc(ancFreq = l, invL = pop_InvL, alpha = pop_alpha))
  pop_p_ZAnc = calcProbZAnc(ZAnc = pop_ZAnc, nPop = length(pop_alpha), logProb = F)
  return(data.frame(ZAnc = pop_ZAnc, p_ZAnc = pop_p_ZAnc, mean_anc = mean_anc,
                    stringsAsFactors = F))
}

# file with list of paths to population ancestry files
pop_files = read.table(file.path("results", prefix, "pop.anc.list"),
                       stringsAsFactors = F, header = F)$V1

# read in population ancestry input files from calc_genomewide_pop_anc_freq.R
pop_anc_list = lapply(pop_files, function(pop) read.table(pop, 
                                                          stringsAsFactors = F))
# combine population ancestry frequencies into a matrix
# where rows are populations and columns are snps
all_anc = t(do.call(cbind, pop_anc_list))

# run K-matrix and ZAnc calculations:
all = make_calcs(pop_anc = all_anc)

# print results
write.table(all, 
            file.path("results", prefix, "pop.ZAnc"),
            sep = "\t", 
            col.names = T, row.names = F, quote = F)


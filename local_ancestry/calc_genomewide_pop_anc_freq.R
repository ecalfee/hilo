#!/usr/bin/env Rscript

# this script takes in a hilo admixed population #
# and a directory where to find the HILOXX.posterior files from ancestry_hmm
# and outputs 3 files into a new 'anc' subdirectory:
# a popN.anc.ind file with ancestry for all individuals in popN individually
# a popN.anc.freq file with mean ancestry for all individuals in popN
# a popN.alpha.ind file with two columns, n (ind HILO ID)
# and alpha (individual genomewide mean ancestry from local ancestry posterior)
# alpha and anc represent the % mexicana ancestry (as opposed to maize)

library(dplyr)

# to run:
# Rscript calc_genomewide_pop_anc_freq.R 366 results/ancestry_hmm/pass2_alloMAIZE/output_noBoot rsults/ancestry_hmm/pass2_alloMAIZE/input/pop366.anc_hmm.ids

# arguments
args = commandArgs(trailingOnly=TRUE)
# population number
N = as.integer(args[1])
# paths to input and output directories
dir_input = args[2]
# individuals file with IDs of samples with local ancestry calls:
ids_file = args[3]
dir_output = paste0(dir_input, "/anc")
if (file.exists(dir_input) & !file.exists(dir_output)){ # make sure output directory exists first (or create)!
  dir.create(file.path(dir_output), recursive = T)
}

# find included individuals from current pop
pop_ids <- read.table(ids_file, stringsAsFactors = F,
                    header = F, sep = "\t")$V1

# function to summarise individual ancestry from posterior output
sample_pos = function(ID, path) {# summarize across posterior of all 3 genotypes 
  # to get estimate of mex. ancestry at each of all SNPs for a HILO individual ID
  as.matrix(read.table(paste0(path, "/", ID, ".posterior"), 
                       header = T, 
                       stringsAsFactors = F)[ , 3:5]) %*% c(0, .5, 1) 
}
sample_pos_nth = function(nth, ID, path) {# take every nth SNP & summarize across posterior of all 3 genotypes 
  # to get estimate of mex. ancestry at each of those SNPs for a HILO individual ID
  as.matrix(read.table(paste0(path, "/", ID, ".posterior"), 
                       header = T, 
                       stringsAsFactors = F)[c(T, rep(F, nth-1)), 3:5]) %*% 
    c(0, .5, 1) 
}
# get ancestry for all individuals
# rows = SNPs; columns = individuals
anc = do.call(cbind, 
                lapply(pop_ids, 
                       function(id) sample_pos(ID = id, path = dir_input)))
# mean mexicana ancestry across all markers per individual
alpha = data.frame(hilo_id = pop_ids, alpha = apply(anc, 2, mean))

# mean population ancestry at each SNP
pop_anc_freq = apply(anc, 1, mean)

# write output files:
# a popN.anc.ind file with ancestry for all individuals in popN individually
write.table(anc, 
            paste0(dir_output, "/pop", N, ".anc.ind"), 
            col.names = F, row.names = F, quote = F, sep = "\t")

# a popN.anc.freq file with mean ancestry for all individuals in popN
write.table(pop_anc_freq, 
            paste0(dir_output, "/pop", N, ".anc.freq"), 
            col.names = F, row.names = F, quote = F, sep = "\t")

# a popN.alpha.ind file with two columns, n (ind HILO ID)
write.table(alpha, 
            paste0(dir_output, "/pop", N, ".alpha.ind"), 
            col.names = T, row.names = F, quote = F, sep = "\t")


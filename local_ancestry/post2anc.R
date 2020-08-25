#!/usr/bin/env Rscript
library(dplyr)

# this script takes in an admixed population ploidy file with HILO IDs
# and a directory where to find the HILOXX.posterior files from ancestry_hmm
# and outputs 3 files into a new 'anc' subdirectory:
# a popN.anc.ind file with ancestry for all individuals in popN individually
# a popN.anc.freq file with mean ancestry for all individuals in popN
# a popN.alpha.ind file with two columns, n (ind HILO ID)
# and alpha (individual genomewide mean ancestry from local ancestry posterior)
# alpha and anc represent the % mexicana ancestry (as opposed to maize)

# load variables from Snakefile
ploidy_file = snakemake@input[["ploidy"]]
dir_input = snakemake@params[["dir_post"]]
output_file_anc_pop = snakemake@output[["pop"]]
output_file_anc_ind = snakemake@output[["ind"]]
output_file_alpha_ind = snakemake@output[["alpha"]]


# find included individuals from current pop
pop_ids <- read.table(ploidy_file, stringsAsFactors = F,
                    header = F, sep = "\t")$V1

# function to summarise individual ancestry from posterior output
sample_pos = function(ID, path) {# summarize across posterior of all 3 genotypes 
  # to get estimate of mex. ancestry at each of all SNPs for a HILO individual ID
  as.matrix(read.table(paste0(path, "/", ID, ".posterior"), 
                       header = T, 
                       stringsAsFactors = F)[ , 3:5]) %*% c(0, .5, 1) 
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
            output_file_anc_ind, 
            col.names = F, row.names = F, quote = F, sep = "\t")

# a popN.anc.freq file with mean ancestry for all individuals in popN
write.table(pop_anc_freq, 
            output_file_anc_pop, 
            col.names = F, row.names = F, quote = F, sep = "\t")

# a popN.alpha.ind file with two columns, n (ind HILO ID)
write.table(alpha, 
            output_file_alpha_ind, 
            col.names = T, row.names = F, quote = F, sep = "\t")


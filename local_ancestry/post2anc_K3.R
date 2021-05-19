#!/usr/bin/env Rscript
library(dplyr)

# this script takes in an admixed population ploidy file with HILO IDs
# and a directory where to find the HILOXX.posterior files from ancestry_hmm
# and outputs 3 files into a new 'anc/ZEA' subdirectory:
# a popN.anc.ind file with ancestry for all individuals in popN individually
# a popN.anc.freq file with mean ancestry for all individuals in popN
# a popN.alpha.ind file with two columns, n (ind HILO ID)
# and alpha (individual genomewide mean ancestry from local ancestry posterior)
# alpha and anc represent the % ancestry (ZEA is mexicana or maize)

# load variables from Snakefile
ploidy_file = snakemake@input[["ploidy"]]
dir_input = snakemake@params[["dir_post"]]
output_file_anc_pop_mex = snakemake@output[["pop_mex"]]
output_file_anc_ind_mex = snakemake@output[["ind_mex"]]
output_file_alpha_ind_mex = snakemake@output[["alpha_mex"]]
output_file_anc_pop_parv = snakemake@output[["pop_parv"]]
output_file_anc_ind_parv = snakemake@output[["ind_parv"]]
output_file_alpha_ind_parv = snakemake@output[["alpha_parv"]]
output_file_anc_pop_maize = snakemake@output[["pop_maize"]]
output_file_anc_ind_maize = snakemake@output[["ind_maize"]]
output_file_alpha_ind_maize = snakemake@output[["alpha_maize"]]



# find included individuals from current pop
pop_ids <- read.table(ploidy_file, stringsAsFactors = F,
                    header = F, sep = "\t")$V1

# genotypes
genotypes = c("X2.0.0", "X1.1.0", "X1.0.1", "X0.2.0", "X0.1.1", "X0.0.2") # used in sample_pos
# assumes order of admixture is mexicana, then parviglumis, then maize:
mexicana_weights <- c(1, .5, .5, 0, 0, 0) # ancestry weighting of columns in order above ('genotypes')
parv_weights <- c(0, .5, 0, 1, .5, 0)
maize_weights <- c(0, 0, .5, 0, .5, 1)

# function to summarise individual ancestry from posterior output
sample_pos = function(ID, path, weights) {# summarize across posterior of all 3 genotypes
  # to get estimate of mex. ancestry at each of all SNPs for a HILO individual ID
  as.matrix(read.table(paste0(path, "/", ID, ".posterior"),
                       header = T,
                       stringsAsFactors = F)[ , genotypes]) %*% weights
}


# get ancestry for all individuals
# in format rows = SNPs; columns = individuals
# from individual posterior files with columns 'genotypes'
anc_mex = do.call(cbind,
                lapply(pop_ids,
                       function(id) sample_pos(ID = id, path = dir_input,
                       weights = mexicana_weights)))
# mean mexicana ancestry across all markers per individual
alpha_mex = data.frame(hilo_id = pop_ids, alpha = apply(anc_mex, 2, mean))

# mean population ancestry at each SNP
pop_anc_freq_mex = apply(anc_mex, 1, mean)

# write output files:
# a popN.anc.ind file with ancestry for all individuals in popN individually
write.table(anc_mex,
            output_file_anc_ind_mex,
            col.names = F, row.names = F, quote = F, sep = "\t")

# a popN.anc.freq file with mean ancestry for all individuals in popN
write.table(pop_anc_freq_mex,
            output_file_anc_pop_mex,
            col.names = F, row.names = F, quote = F, sep = "\t")

# a popN.alpha.ind file with two columns, n (ind HILO ID)
write.table(alpha_mex,
            output_file_alpha_ind_mex,
            col.names = T, row.names = F, quote = F, sep = "\t")


# now repeat for maize ancestry estimates:
anc_maize = do.call(cbind,
                lapply(pop_ids,
                       function(id) sample_pos(ID = id, path = dir_input,
                       weights = maize_weights)))
# mean maize ancestry across all markers per individual
alpha_maize = data.frame(hilo_id = pop_ids, alpha = apply(anc_maize, 2, mean))

# mean population ancestry at each SNP
pop_anc_freq_maize = apply(anc_maize, 1, mean)

# write output files:
# a popN.anc.ind file with ancestry for all individuals in popN individually
write.table(anc_maize,
            output_file_anc_ind_maize,
            col.names = F, row.names = F, quote = F, sep = "\t")

# a popN.anc.freq file with mean ancestry for all individuals in popN
write.table(pop_anc_freq_maize,
            output_file_anc_pop_maize,
            col.names = F, row.names = F, quote = F, sep = "\t")

# a popN.alpha.ind file with two columns, n (ind HILO ID)
write.table(alpha_maize,
            output_file_alpha_ind_maize,
            col.names = T, row.names = F, quote = F, sep = "\t")

# finally, repeat for parviglumis

anc_parv = do.call(cbind,
                lapply(pop_ids,
                       function(id) sample_pos(ID = id, path = dir_input,
                       weights = parv_weights)))
# mean parv ancestry across all markers per individual
alpha_parv = data.frame(hilo_id = pop_ids, alpha = apply(anc_parv, 2, mean))

# mean population ancestry at each SNP
pop_anc_freq_parv = apply(anc_parv, 1, mean)

# write output files:
# a popN.anc.ind file with ancestry for all individuals in popN individually
write.table(anc_parv,
            output_file_anc_ind_parv,
            col.names = F, row.names = F, quote = F, sep = "\t")

# a popN.anc.freq file with mean ancestry for all individuals in popN
write.table(pop_anc_freq_parv,
            output_file_anc_pop_parv,
            col.names = F, row.names = F, quote = F, sep = "\t")

# a popN.alpha.ind file with two columns, n (ind HILO ID)
write.table(alpha_parv,
            output_file_alpha_ind_parv,
            col.names = T, row.names = F, quote = F, sep = "\t")

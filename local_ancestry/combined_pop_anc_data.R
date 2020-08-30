#!/usr/bin/env Rscript
library(dplyr)

# this script takes in population ancestry frequencies for all sympatric pops 
# for maize or mexicana and saves 
# (1) a summary R dataframe (SNPs X Pops) and 
# (2) a combined sample allele frequency file for all maize/mexicana

# load variables from Snakefile
meta_file = snakemake@input[[1]]
# meta_file = "samples/HILO_MAIZE55_meta.RData"
k2_file = snakemake@input[[2]]
dir_anc = snakemake@params[["dir_anc"]]
dir_ploidy = snakemake@params[["dir_ploidy"]]
zea = snakemake@parmas[["zea"]]
output_pops = snakemake@output[["pops"]]
output_combined = snakemake@output[["combined"]]
output_meta = snakemake@output[["meta"]]

# load metadata
load(meta_file) # meta

k2 = read.table(k2_file, sep ="\t", header = T, stringsAsFactors = F)

pops = unique(filter(meta, symp_allo == "sympatric" & zea == zea)$popN) 

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


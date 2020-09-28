#!/usr/bin/env Rscript
library(dplyr)

# this script takes in an admixed population ploidy file with HILO IDs
# and a directory where to find the HILOXX.posterior files from ancestry_hmm
# and outputs 2 files into the new 'MAP' subdirectory:
# a MAP/popN.anc.ind file with the highest posterior probability ancestry state for that ind. at each site along the genome
# 1 = homozygous mexicana, 0.5 = het, 0 = homozygous maize
# a MAP/popN.prob.ind file with the posterior probability for that ancestry state

# load variables from Snakefile
ploidy_file = snakemake@input[["ploidy"]]
dir_input = snakemake@params[["dir_post"]]
output_file_map_anc = snakemake@output[["anc"]]
output_file_map_prob = snakemake@output[["prob"]]


# find included individuals from current pop
pop_ids <- read.table(ploidy_file, stringsAsFactors = F,
                    header = F, sep = "\t")$V1
# function to read in posterior probabilities from ancestry_hmm output
# column order: homozygous maize, heterozygous, homozygous mexicana ancestry
get_post = function(ID, path) {
  post = read.table(paste0(path, "/", ID, ".posterior"), 
                    header = T, 
                    stringsAsFactors = F)[ , 3:5]
}
# function to get MAP individual ancestry and probability from posterior output
# note: for true ties, MAP is set randomly to one of the tied states
get_map = function(post){
  prob = apply(post, 1, max)
  anc_values = c(0,0.5,1)
  max_anc = apply(post, 1, function(row)
    anc_values[sample(which(row == max(row)), 1)]) # sampling 1 value picks randomly in case of ties
  return(data.frame(anc = max_anc, prob = prob, stringsAsFactors = F))
}

# test function get_map()
# post1 <- data.frame(AA = c(1,.9,0.2,0.1,0.5), 
#                     Aa = c(0,.05,0.4,0.4,0.5), 
#                     aa = c(0,0.05,0.4,0.5,0))
# get_map(post1)

# get MAP ancestry for all individuals
# rows = SNPs; columns = individuals
map_all = lapply(pop_ids, 
             function(id) get_map(get_post(ID = id, path = dir_input)))

# MAP ancestry state
anc = do.call(cbind, lapply(map_all, function(x) x$anc))

# MAP ancestry probability
prob = do.call(cbind, lapply(map_all, function(x) x$prob))

# write output files:
write.table(anc, 
            output_file_map_anc, 
            col.names = F, row.names = F, quote = F, sep = "\t")
write.table(prob, 
            output_file_map_prob, 
            col.names = F, row.names = F, quote = F, sep = "\t")
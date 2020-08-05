#!/usr/bin/env Rscript

# This script assigns the labels 'maize' 'mexicana' 'parviglumis'
# to the ancestries that come out of NGSadmix, based
# on which allopatric sample (zea subspecies) has the most of that ancestry
library(dplyr)
library(tidyr)

# load variables from Snakefile
K = snakemake@params["k"]
# K = 2
admix_file = snakemake@input["admix"]
file_out = snakemake@output["anc"]
load(snakemake@input[["meta"]]) # sample metadata - samples are in the same order as NGSAdmix results
# load("samples/HILO_MAIZE55_meta.RData")


ancestries <- c("maize", "mexicana", "parviglumis")[1:K]

admix <- read.table(admix_file)
colnames(admix) <- paste0("anc", 1:K) #c("anc1", "anc2")

# join bams and admix by position (CAUTION - bam list order and admix results MUST MATCH!)
d <- bind_cols(meta, admix)  %>%
  arrange(., popN) %>%
  arrange(., zea) %>%
  arrange(., symp_allo)


which_anc <- data.frame(ancestry = colnames(admix),
                   ancestry_label = sapply(colnames(admix), function(x) # for each of K ancestries,
                       names(which.max(tapply(d$x[d$symp_allo == "allopatric"], 
                                              d$zea[d$symp_allo == "allopatric"], 
                                              mean)))), # which zea allopatric population, on avg, has the highest % of that ancestry?
                   stringsAsFactors = F)

d %>% 
  tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  left_join(., which_anc, by = "ancestry") %>%
  dplyr::select(-ancestry) %>%
  tidyr::spread(., ancestry_label, p) %>%
  dplyr::select(ID, ancestries) %>%
  write.table(., file_out,
              quote = F,
              sep = "\t", col.names = T, row.names = F)
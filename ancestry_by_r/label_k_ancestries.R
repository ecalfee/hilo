#!/usr/bin/env Rscript

# This script assigns the labels 'maize' 'mexicana' 'parviglumis'
# to the ancestries that come out of NGSadmix, based
# on which allopatric sample (zea subspecies) has the most of that ancestry
library(dplyr)
library(tidyr)

# load variables from Snakefile
K = snakemake@params[["k"]]
# K = 3
admix_file = snakemake@input[["admix"]]
# admix_file = "ancestry_by_r/results/bootstrap_1cM/HILO_MAIZE55_PARV50/cd5_5/K3/boot60.qopt"
file_out = snakemake@output[["anc"]]
# file_out = "ancestry_by_r/results/bootstrap_1cM/HILO_MAIZE55_PARV50/cd5_5/K3/boot60.anc"
load(snakemake@input[["meta"]]) # sample metadata - samples are in the same order as NGSAdmix results
# load("samples/HILO_MAIZE55_PARV50_meta.RData")


ancestries <- c("maize", "mexicana", "parviglumis")[1:K]

admix <- read.table(admix_file)
colnames(admix) <- paste0("anc", 1:K) #c("anc1", "anc2")

# join bams and admix by position (CAUTION - bam list order and admix results MUST MATCH!)
d <- bind_cols(meta, admix)  %>%
  arrange(., popN) %>%
  arrange(., zea) %>%
  arrange(., symp_allo) %>%
  tidyr::pivot_longer(data = ., cols = colnames(admix), # separate out the K ancestries to each have their own row
                      names_to = "ancestry", values_to = "p")

# for each of k unlabelled ancestries, find the zea subspecies that is the best match
which_anc <- d %>%
  dplyr::mutate(symp_allo = ifelse(zea == "parviglumis", "allopatric", symp_allo)) %>%
  filter(., symp_allo == "allopatric") %>% # only include allopatric samples
  group_by(zea, ancestry) %>%  
  summarise(p = mean(p)) %>% # for each zea subspecies, what is the mean proportion p of each ancestry?
  group_by(ancestry) %>%
  summarise(ancestry_label = zea[which.max(p)]) # which subspecies has the most of each ancestry? assign that zea as the ancestry label

if (length(unique(which_anc$ancestry_label)) == K){ # ancestries unambiguously assigned
  d_out = d %>%
    left_join(., which_anc, by = "ancestry") %>%
    dplyr::select(-ancestry) %>%
    tidyr::spread(., ancestry_label, p) %>%
    dplyr::select(ID, ancestries)
}else {# ancestries can't be unambiguously assigned
  print("error: ancestries could not be unambiguously assigned")
  which_anc
  d_out = dplyr::select(d, ID)
  for (a in ancestries){
    d_out[ , a] = NA # set all ancestry values to NA
  }
}

# write output
write.table(d_out, file_out,
              quote = F,
              sep = "\t", col.names = T, row.names = F)

#!/usr/bin/env Rscript
# working directory is hilo/
# summarising NGSAdmix global ancestry results 
# for low-coverage maize/mexicana individuals 
# in HILO adaptation project
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)

# load variables from Snakefile
K = snakemake@params[["k"]]
# K = 2
admix_file = as.numeric(snakemake@input[["admix"]])
# admix_file = "global_ancestry/results/NGSAdmix/HILO_MAIZE55/K2.qopt"
alphas_out = snakemake@output[["alphas"]]
# alphas_out = "global_ancestry/results/NGSAdmix/HILO_MAIZE55/K2_alphas_by_symp_pop.txt"
ind_out = snakemake@output[["ind"]]
# ind_out = "global_ancestry/results/NGSAdmix/HILO_MAIZE55/K2_alphas_by_ind.RData"
load(snakemake@input[["meta"]]) # sample metadata, including estimated sequence coverage
# load("samples/HILO_MAIZE55_meta.RData")
# note: # metadata in same order as samples in admixture

# NGSAdmix data
#ancestries <- c("maize", "mexicana", "parviglumis")[1:K]
admix <- read.table(admix_file)
colnames(admix) <- paste0("anc", 1:K) #c("anc1", "anc2")

# join bams and admix by position (CAUTION - bam list order and admix results MUST MATCH!)
d <- bind_cols(meta, admix)  %>%
  arrange(., popN) %>%
  arrange(., zea) %>%
  arrange(., symp_allo) %>%
  arrange(., group, ELEVATION)

# which unlabelled ancestry maps onto maize and which onto mexicana?
anc <- d %>%
  pivot_longer(., cols = colnames(admix), names_to = "K_anc", values_to = "p_anc") %>%
  group_by(zea) %>%
  summarise(anc_label = K_anc[which.max(p_anc)])

# label NGSAdmix ancestry proportions as alpha_mex and alpha_maize
d_admix2 <- rename(d, mexicana = anc$anc_label[anc$zea == "mexicana"], 
                  maize = anc$anc_label[anc$zea == "maize"])


# calculate mean proportion mexicana ancestry per sympatric population for ancestry_hmm input
alphasByPop = d_admix2 %>%
  #filter(., symp_allo == "sympatric" & est_coverage >= min_coverage) %>%
  filter(., symp_allo == "sympatric") %>%
  group_by(., popN) %>%
  summarise(., alpha_maize = mean(maize), 
            alpha_mex = mean(mexicana),
            n_global_ancestry = n(),
            n_local_ancestry = sum(est_coverage >= 0.5)) # number included individuals per group

# write file with mean population ancestry (for prior in ancestry_hmm)
write.table(format(alphasByPop, scientific = F),  # write table so I can pull population global ancestry values for ancestry_hmm
            file.path(alphas_out),
            col.names = T, row.names = F, quote = F, sep = "\t")

# also save estimated admixture proportions for individuals
save(list = c("d_admix2"), file = ind_out)

# # lets look at allele frequency estimates for each ancestry
# #freqs <- read.table("../data/NGSadmix/pass1/K2_pruned_all.fopt.gz") 
# freqs <- read.table(paste0(file_prefix, ".fopt.gz"))
# colnames(freqs) <- paste0("anc", 1:K)
# #freqs$diff <- abs(freqs$anc1 - freqs$anc2)
# mean_allele_freq <- apply(freqs, 1, mean)
# hetSubpops <- sapply(1:K, function(i) 
#   mean(2 * freqs[ , paste0("anc", i)] * 
#          (1 - freqs[ , paste0("anc", i)])))
# #het1 <- mean(freqs$anc1 * (1-freqs$anc1))
# #het2 <- mean(freqs$anc2 * (1-freqs$anc2))
# # a better total approximate would be to use the MAF of the entire sample estimated in ANGSD
# # rather than the mean of equal sized hypothetical pop 'ancestries'
# hetTot <- mean(2 * mean_allele_freq * (1 - mean_allele_freq))
# # differentiation is fairly low between the two identified genetic clusters
# Fst <- (hetTot - mean(hetSubpops))/hetTot
# Fst

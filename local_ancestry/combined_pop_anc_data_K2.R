#!/usr/bin/env Rscript
library(dplyr)

# this script takes in population ancestry frequencies for all sympatric pops
# for maize or mexicana and saves
# (1) a summary R dataframe (SNPs X Pops) and
# (2) a combined sample allele frequency file for all maize/mexicana

# load variables from Snakefile
meta_file = snakemake@input[[1]]
# meta_file = "samples/HILO_MAIZE55_meta.RData"
tracts_file = snakemake@input[[2]]
# tracts_file = "local_ancestry/results/thinnedSNPs/HILO_MAIZE55/whole_genome.bed"
dir_anc = snakemake@params[["dir_anc"]]
# dir_anc = "local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc"
dir_ploidy = snakemake@params[["dir_ploidy"]]
# dir_ploidy = "local_ancestry/results/ancestry_hmm/HILO_MAIZE55/input"
ZEA = snakemake@params[["zea"]]
# ZEA = "mexicana"
output_pop_anc = snakemake@output[["pop_anc"]]
output_meta_pop = snakemake@output[["meta_pop"]]
output_mexicana_anc = snakemake@output[["mexicana_anc"]]
output_maize_anc = snakemake@output[["maize_anc"]]

# load metadata
load(meta_file) # meta

tracts = read.table(tracts_file, sep = "\t", header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("chr", "start", "end", "pos"))

# which pops are included?
meta_pops = dplyr::filter(meta, symp_allo == "sympatric" & zea == ZEA) %>%
  dplyr::select(popN, zea, symp_allo, group, GEOCTY, LOCALITY, ELEVATION, LAT, LON) %>%
  arrange(ELEVATION) %>%
  distinct() %>%
  mutate(pop = paste0("pop", popN),
         n_local_ancestry = sapply(popN, function(i)
           nrow(read.table(paste0(dir_ploidy, "/pop", i, ".ploidy"),
                             stringsAsFactors = F,
                             header = F, sep = "\t"))))

# get ancestry for all individuals
# rows = SNPs; columns = individuals
# mexicana ancestry
anc_mexicana = do.call(cbind,
                lapply(meta_pops$popN,
                       function(i) read.table(paste0(dir_anc, "/mexicana/pop", i, ".anc.freq"),
                                          header = F, stringsAsFactors = F)$V1))
colnames(anc_mexicana) = meta_pops$pop
# maize ancestry
anc_maize = do.call(cbind,
                lapply(meta_pops$popN,
                       function(i) read.table(paste0(dir_anc, "/maize/pop", i, ".anc.freq"),
                                          header = F, stringsAsFactors = F)$V1))
colnames(anc_maize) = meta_pops$pop

# mean mexicana/maize ancestry across all markers per pop
meta_pops$alpha_local_ancestry_mexicana = apply(anc_mexicana, 2, mean)
meta_pops$alpha_local_ancestry_maize = apply(anc_maize, 2, mean)

# mean ancestry at each SNP, combined across all indidivuals
anc_mexicana_mean = (anc_mexicana %*% meta_pops$n_local_ancestry)/sum(meta_pops$n_local_ancestry)
anc_maize_mean = (anc_maize %*% meta_pops$n_local_ancestry)/sum(meta_pops$n_local_ancestry)

anc = list()
anc[["mexicana"]] = anc_mexicana
anc[["maize"]] = anc_maize
anc_mean = list()
anc_mean[["mexicana"]] = anc_mexicana_mean
anc_mean[["maize"]] = anc_maize_mean

# write output files:

# RData files
save(list = "meta_pops", file = output_meta_pop)
save(list = c("anc", "anc_mean"), file = output_pop_anc)

# write bedfiles -- 1 per ancestry -- with population ancestry frequencies for each tract
write.table(tracts %>%
                mutate(., anc_freq = anc_mexicana_mean),
            file = output_mexicana_anc,
            col.names = T, row.names = F, quote = F, sep = "\t")
write.table(tracts %>%
                mutate(., anc_freq = anc_maize_mean),
                file = output_maize_anc,
                col.names = T, row.names = F, quote = F, sep = "\t")

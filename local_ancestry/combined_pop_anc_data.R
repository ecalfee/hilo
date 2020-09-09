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
ZEA = snakemake@parmas[["zea"]]
# ZEA = "mexicana"
output_pop_anc = snakemake@output[["pop_anc"]]
output_combined_anc = snakemake@output[["combined_anc"]]
output_meta_pop = snakemake@output[["meta_pop"]]

# load metadata
load(meta_file) # meta

tracts = read.table(tracts_file, sep = "\t", header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("chr", "start", "end", "pos"))

# which pops are included?
meta_pops = dplyr::filter(meta, symp_allo == "sympatric" & zea == ZEA) %>%
  dplyr::select(popN, zea, symp_allo, group, GEOCTY, LOCALITY, ELEVATION, LAT, LON) %>%
  arrange(ELEVATION) %>%
  distinct() %>%
  mutate(pop = paste0("pop", popN)) %>%
  mutate(n_local_ancestry = nrow(read.table(paste0(dir_ploidy, "/pop", popN, ".ploidy"), 
                             stringsAsFactors = F,
                             header = F, sep = "\t")))

# get ancestry for all individuals
# rows = SNPs; columns = individuals
anc = do.call(cbind, 
                lapply(meta_pops$popN, 
                       function(i) read.table(paste0(dir_anc, "/pop", i, ".anc.freq"),
                                          header = F, stringsAsFactors = F)$V1))
colnames(anc) = meta_pops$pop

# mean mexicana ancestry across all markers per pop
meta_pops$alpha_local_ancestry = apply(anc, 2, mean)

# mean ancestry at each SNP, combined across all indidivuals
combined_anc_freq = (anc %*% meta_pops$n_local_ancestry)/sum(meta_pops$n_local_ancestry)
tracts$anc_freq = combined_anc_freq

# write output files:
save(list = "meta_pops", file = output_meta_pop)
save(list = "anc", file = output_pop_anc)
write.table(tracts, # bed file with combined ancestry frequencies
            file = output_combined_anc, 
            col.names = T, row.names = F, quote = F, sep = "\t")
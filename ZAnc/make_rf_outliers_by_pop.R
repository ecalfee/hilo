#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(bedr)


# this script identifies high introgression ancestry outlier regions
# in an individual population and shared across pops
# and outputs regions files for outliers
# from focal pop only (1pop)
# and from focal pop + at least 3 other pops (4pop)
# to be used with angsd -rf

# load variables from Snakefile
bed_sites = snakemake@input[["bed_sites"]]
# bed_sites = "local_ancestry/results/thinnedSNPs/HILO_MAIZE55/whole_genome.bed"
zea = snakemake@params[["zea"]]
# zea = "maize"
anc_file = snakemake@input[["anc"]]
# anc_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pops.anc.RData")
meta_file = snakemake@input[["meta"]]
# meta_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pop.meta.RData")
dir_out = snakemake@params[["dir_out"]]
# dir_out = "ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot"
focal_pop = snakemake@params[["focal_pop"]]
# focal_pop = "pop362"

# load data
load(meta_file)
load(anc_file)

sites <- read.table(bed_sites, header = F, stringsAsFactors = F, sep = "\t") %>%
  data.table::setnames(c("chr", "start", "end", "length"))

if (zea == "mexicana"){
    anc = 1 - anc # always count minor (introgressed) ancestry
    meta_pops$alpha_local_ancestry = 1 - meta_pops$alpha_local_ancestry
}
  
# what is 2sd above the mean introgressed ancestry threshold for each pop?
meta_pops$sd2 <- apply(anc, 2, mean) + 2*apply(anc, 2, sd)
  
# find pop outliers in observed data
anc_outliers <- data.frame(anc, stringsAsFactors = F) %>%
    cbind(., sites) %>%
    tidyr::pivot_longer(., cols = colnames(anc), names_to = "pop", values_to = "anc") %>%
    left_join(., meta_pops, by = "pop") %>%
    mutate(top_sd2 = anc > sd2) %>%
    arrange(ELEVATION)
  
# add a column for how many populations are an outlier at that position
# and filter to only include outliers that involve the focal population
anc_outliers_by_pop <- anc_outliers %>%
    group_by(chr, start, end) %>%
    summarise(outlier_pops = sum(top_sd2)) %>%
    ungroup() %>%
  left_join(anc_outliers, ., by = c("chr", "start", "end")) %>%
  filter(pop == focal_pop & top_sd2) %>% # must be an outlier in focal pop for that row
  mutate(region = paste0(chr, ":", start + 1, "-", end)) # format for angsd regions file

# print regions files (rf) for single pop outlier regions
anc_outliers_by_pop %>%
    filter(outlier_pops == 1) %>%
    dplyr::select(region) %>% # only print region
    write.table(., file = paste0(dir_out, "/", focal_pop, ".1pop.outliers.regions"),
                col.names = F, row.names = F, sep = "\t", quote = F)

# print regions files (rf) for 4+ pop outlier regions
anc_outliers_by_pop %>%
  filter(outlier_pops >= 4) %>%
  dplyr::select(region) %>% # only print region
  write.table(., file = paste0(dir_out, "/", focal_pop, ".4pop.outliers.regions"),
              col.names = F, row.names = F, sep = "\t", quote = F)


# print bed file for single pop outlier regions
anc_outliers_by_pop %>%
  filter(outlier_pops == 1) %>%
  dplyr::select(chr, start, end) %>% # only print region
  write.table(., file = paste0(dir_out, "/", focal_pop, ".1pop.outliers.bed"),
              col.names = F, row.names = F, sep = "\t", quote = F)

# print bed file for 4+ pop outlier regions
anc_outliers_by_pop %>%
  filter(outlier_pops >= 4) %>%
  dplyr::select(chr, start, end) %>% # only print region
  write.table(., file = paste0(dir_out, "/", focal_pop, ".4pop.outliers.bed"),
              col.names = F, row.names = F, sep = "\t", quote = F)

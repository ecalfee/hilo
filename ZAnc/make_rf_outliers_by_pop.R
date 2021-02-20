#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)

# this script identifies high introgression ancestry outlier regions
# in an individual population and shared across pops
# and outputs regions files for outliers
# from focal pop only (1pop)
# and from focal pop + at least 3 other pops (4pop)
# to be used with angsd -rf

# load variables from Snakefile
bed_sites = snakemake@input[["bed_sites"]]
# bed_sites = "local_ancestry/results/thinnedSNPs/HILO_MAIZE55/whole_genome.bed"
anc_maize = snakemake@input[["anc_maize"]]
# anc_maize = "local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/maize.pops.anc.RData"
meta_maize = snakemake@input[["meta_maize"]]
# meta_maize = "local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/maize.pop.meta.RData"
anc_mexicana = snakemake@input[["anc_mexicana"]]
# anc_mexicana = "local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/mexicana.pops.anc.RData"
meta_mexicana = snakemake@input[["meta_mexicana"]]
# meta_mexicana = "local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/mexicana.pop.meta.RData"
dir_out = snakemake@params[["dir_out"]]
# dir_out = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot") 
focal_pop = snakemake@params[["focal_pop"]]
# focal_pop = "pop362"
meta_file = snakemake@input[["meta"]]
# meta_file = "samples/HILO_MAIZE55_meta.RData"

# is the focal population sympatric maize or mexicana?
load(meta_file)
meta_sympatric = meta %>%
  filter(symp_allo == "sympatric") %>%
  dplyr::select(popN, zea, symp_allo, LOCALITY, ELEVATION) %>%
  filter(!duplicated(.)) %>%
  mutate(pop = paste0("pop", popN))

zea = meta_sympatric$zea[meta_sympatric$pop == focal_pop]

# load ancestry and population metadata files
# based on whether the focal population is sympatric
# maize or mexicana
if (zea == "maize"){
  load(anc_maize)
  load(meta_maize)
}

if (zea == "mexicana"){
  load(anc_mexicana)
  load(meta_mexicana)
  anc = 1 - anc # always count minor (introgressed) ancestry
  meta_pops$alpha_local_ancestry = 1 - meta_pops$alpha_local_ancestry
}

sites <- read.table(bed_sites, header = F, stringsAsFactors = F, sep = "\t") %>%
  data.table::setnames(c("chr", "start", "end", "length"))

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

#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
# this script plots timing of admixture estimated from ancestry_hmm for different Nes
# and correlations between local ancestry estimates for different Nes

# load variables from Snakefile
prefix = snakemake@params[["prefix"]]
# prefix = "HILO_MAIZE55_PARV50"

png_times = snakemake@output[["png_times"]]
# png_times = paste0("local_ancestry/plots/", prefix, "_K3_sensitivity_to_Ne_admix_times.png")
png_times_lzw = snakemake@output[["png_times_lzw"]]
# png_times_lzw = paste0("../hilo_manuscript/figures_supp/", prefix, "_K3_sensitivity_to_Ne_admix_times.tif")
png_anc = snakemake@output[["png_anc"]]
# png_anc = paste0("local_ancestry/plots/", prefix, "_K3_sensitivity_to_Ne_local_ancestry.png")
png_anc_lzw = snakemake@output[["png_anc_lzw"]]
# png_anc_lzw = paste0("../hilo_manuscript/figures_supp/", prefix, "_K3_sensitivity_to_Ne_local_ancestry.tif")

# variables
Nes = c(1000, 10000, 100000)
zea = c("maize", "mexicana")
ancestries = c("maize", "mexicana", "parv")
YESNO = "no"

# load data
source(snakemake@params[["colors"]]) # plotting colors
# source("colors.R")
load(paste0("samples/", prefix, "_meta.RData"))

times = list()

for (z in zea){
  # What to do with Nes???
  # load meta data
  load(paste0("local_ancestry/results/ancestry_hmm/",
              prefix, "/K3/Ne", Ne, "_", YESNO,
              "Boot/anc/", z, ".pop.meta.RData"))
  times[[z]] = do.call(bind_cols, lapply(meta_pops$pop, function(p)
    read.table(paste0("local_ancestry/results/ancestry_hmm/",
                      prefix, "/K3/Ne", Ne, "_", YESNO,
                      "Boot/", p, ".times"),
               header = T, stringsAsFactors = F, sep = " ") %>%
      dplyr::mutate(pop = p))) %>%
    left_join(., dplyr::select(meta_pops, pop, zea, LOCALITY, ELEVATION), by = "pop")
}




# p_times
ggsave(file = png_times,
       plot = p_times,
       device = "png",
       width = 5, height = 4,
       units = "in", dpi = 300)

ggsave(file = png_times_lzw,
       plot = p_times,
       device = "tiff",
       width = 5, height = 4,
       compression = "lzw", type = "cairo",
       units = "in", dpi = 300)

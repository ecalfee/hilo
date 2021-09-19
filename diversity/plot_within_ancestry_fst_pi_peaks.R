#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

# script to summarise diversity (pi) within populations
# and pairwise differentiation (fst) between local maize-mexicana populations
# for homozygous mexicana ancestry windows within high introgression outlier peaks
# in sympatric maize

# get snakemake variables
source(snakemake@input[["colors"]])
# source("colors.R")
fst_peaks_file = snakemake@input[["fst_peaks"]]
# fst_peaks_file = "diversity/results/fst/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/HOMOZYG/summary_pop_pairs_fst.mexicana_ancestry.peaks.allChr.txt"
fst_genomewide_file = snakemake@input[["fst_genomewide"]]
# fst_genomewide_file = "diversity/results/fst/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/HOMOZYG/summary_pop_pairs_fst.allChr.txt"
pi_peaks_file = snakemake@input[["pi_peaks"]]
# pi_peaks_file = "diversity/results/pi/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/HOMOZYG/summary_pop_pi.mexicana_ancestry.peaks.allChr.txt"
pi_genomewide_file = snakemake@input[["pi_genomewide"]]
# pi_genomewide_file = "diversity/results/pi/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/HOMOZYG/summary_pop_pi.allChr.txt"
load(snakemake@input[["meta"]])
# load("samples/HILO_MAIZE55_PARV50_meta.RData")

# output plots
png_fst_mexicana_anc = snakemake@output[["png_fst_mex"]]
# png_fst_mexicana_anc = "diversity/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_fst_within_mexicana_ancestry_peaks.png"
png_fst_mexicana_anc_lzw = snakemake@output[["png_fst_mex_lzw"]]
# png_fst_mexicana_anc_lzw = "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_fst_within_mexicana_ancestry_peaks.tif"

png_pi_anc = snakemake@output[["png_pi_anc"]]
# png_pi_anc = "diversity/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_pi_within_ancestry.png"
png_pi_anc_lzw = snakemake@output[["png_pi_anc_lzw"]]
# png_pi_anc_lzw = "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_pi_within_ancestry.tif"

meta_pops = meta %>%
  dplyr::select(popN, zea, symp_allo, group, LOCALITY, ELEVATION, LAT, LON) %>%
  dplyr::distinct() %>%
  dplyr::mutate(., pop = paste0("pop", popN)) %>%
  arrange(zea, ELEVATION) %>%
  filter(symp_allo == "sympatric")

fst_genomewide <- read.table(fst_genomewide_file, sep = " ", header = F) %>%
  data.table::setnames(c("file_name", "fst")) %>%
  tidyr::separate(data = .,
                  col = file_name, 
                  into = c("v1", "v2", "v3", "v4", "v5", "v6", "v7", "ancestry", "pops"), 
                  sep = "/") %>%
  tidyr::separate(data = .,
                  col = pops,
                  into = c("pop1", "pop2", "v8", "v9", "v10"),
                  sep = "[.]") %>%
  dplyr::mutate(genomic_region = "genomewide") %>%
  dplyr::select(ancestry, pop1, pop2, genomic_region, fst)

fst_peaks <- read.table(fst_peaks_file, sep = " ", header = F) %>%
  data.table::setnames(c("file_name", "fst")) %>%
  tidyr::separate(data = .,
                  col = file_name, 
                  into = c("v1", "v2", "v3", "v4", "v5", "v6", "v7", "ancestry", "pops"), 
                  sep = "/") %>%
  tidyr::separate(data = .,
                  col = pops,
                  into = c("pop1", "pop2", "genomic_region", "v8", "v9", "v10", "v11"),
                  sep = "[.]") %>%
  dplyr::select(ancestry, pop1, pop2, genomic_region, fst)

# add in metadata for both populations being compared,
# and reverse comparison so a-b and b-a are both in dataset
fst <- bind_rows(dplyr::mutate(fst_genomewide, comparison_order = 1), # because fst is symmetric, add in all 'reverse' comparisons
            dplyr::mutate(fst_genomewide, comparison_order = 2) %>%
            dplyr::rename(., pop1 = pop2, pop2 = pop1),
            fst_peaks) %>%
  left_join(., 
            meta_pops, 
            by = c("pop1"="pop")) %>%
  left_join(., 
            meta_pops, 
            by = c("pop2" = "pop"),
            suffix = c("", ".pop2")) %>%
  dplyr::rename(popN.pop1 = popN, 
                zea.pop1 = zea,
                symp_allo.pop1 = symp_allo,
                group.pop1 = group,
                LOCALITY.pop1 = LOCALITY,
                ELEVATION.pop1 = ELEVATION,
                LAT.pop1 = LAT,
                LON.pop1 = LON)

pi_genomewide <- read.table(pi_genomewide_file, sep = " ", header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("file_name", "chr", "pi", "nSites")) %>%
  dplyr::filter(chr != "Chr") %>% # get rid of headers
  dplyr::mutate(pi = as.numeric(pi),
                nSites= as.numeric(nSites)) %>%
  group_by(file_name) %>%
  summarise(pi = sum(pi)/sum(nSites)) %>% # get mean pi across all chromosomes
  tidyr::separate(data = .,
                  col = file_name, 
                  into = c("v1", "v2", "v3", "v4", "v5", "v6", "v7", "ancestry", "label"), 
                  sep = "/") %>%
  tidyr::separate(data = .,
                  col = label,
                  into = c("pop", "v8", "v9", "v10"),
                  sep = "[.]") %>%
  dplyr::mutate(genomic_region = "genomewide") %>%
  dplyr::select(ancestry, pop, genomic_region, pi)

pi_peaks <- read.table(pi_peaks_file, sep = " ", header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("file_name", "chr", "pi", "nSites")) %>%
  dplyr::filter(chr != "Chr") %>% # get rid of headers
  dplyr::mutate(pi = as.numeric(pi),
                nSites= as.numeric(nSites)) %>%
  group_by(file_name) %>%
  summarise(pi = sum(pi)/sum(nSites)) %>% # get mean pi across all chromosomes
  tidyr::separate(data = .,
                  col = file_name, 
                  into = c("v1", "v2", "v3", "v4", "v5", "v6", "v7", "ancestry", "label"), 
                  sep = "/") %>%
  tidyr::separate(data = .,
                  col = label,
                  into = c("pop", "genomic_region", "outliers_pop", "v8", "v9", "v10"),
                  sep = "[.]") %>%
  dplyr::select(ancestry, pop, genomic_region, pi)

pi = bind_rows(pi_genomewide, pi_peaks) %>%
  left_join(., meta_pops, by = "pop")


# ---------- plot data ---------- #
shapes_peaks <- c(17, 16, 0)
names(shapes_peaks) <- c("1pop", "4pop", "genomewide")
labels_peaks <- c("1 population\nintrogression peaks", "4+ population\nintrogression peaks", "all ancestry tracts\ngenomewide")


p_fst_mex <- fst %>%
  filter(LOCALITY.pop1 == LOCALITY.pop2 & 
           ancestry == "mexicana" &
           zea.pop1 == "maize") %>%
  mutate(LOCALITY.pop1 = reorder(LOCALITY.pop1, ELEVATION.pop1)) %>%
  ggplot(aes(x = LOCALITY.pop1, shape = genomic_region, y = fst)) +
  geom_point(color = col_maize_mex_parv[["maize"]]) +
  coord_flip() +
  theme_light() +
  scale_shape_manual(values = shapes_peaks, 
                     labels = labels_peaks_maize) +
  labs(#title = "Fst within mexicana ancestry\nbetween local population pairs", 
       x = "Location", 
       y = expression(F[ST]), 
       color = NULL, 
       shape = NULL) +
  theme(legend.key.size = unit(10, "mm"))
#p_fst_mex
ggsave(filename = png_fst_mexicana_anc,
       plot = p_fst_mex,
       height = 3, width = 5.5, units = "in", 
       device = "png", dpi = 300)

ggsave(filename = png_fst_mexicana_anc_lzw,
       plot = p_fst_mex,
       height = 3, width = 5.5, units = "in", 
       device = "tiff", dpi = 300,
       compression = "lzw", type = "cairo")

# plot diversity (pi) for maize/mexicana ancestry tracts
# found in sympatric maize and mexicana populations
p_pi_both <- pi %>%
  mutate(LOCALITY = reorder(LOCALITY, ELEVATION)) %>%  
  filter(ancestry %in% c("maize", "mexicana")) %>%
  mutate(ancestry = paste(ancestry, "ancestry")) %>%
  ggplot(aes(x = LOCALITY, color = zea, shape = genomic_region, y = pi)) +
  geom_point() +
  coord_flip() +
  theme_light() +
  facet_grid(ancestry ~ .) +
  scale_shape_manual(values = shapes_peaks,
                     labels = labels_peaks) +
  scale_color_manual(values = col_maize_mex_parv,
                     labels = paste("sympatric", names(col_maize_mex_parv))) +
  labs(
    x = "Location of sampled population", 
    y = expression(pi), 
    color = "Sample", 
    shape = "Ancestry tracts") +
  theme(legend.key.size = unit(10, "mm"))
p_pi_both

ggsave(filename = png_pi_anc,
       plot = p_pi_both,
       height = 6, width = 5.5, units = "in", 
       device = "png", dpi = 300)

ggsave(filename = png_pi_anc_lzw,
       plot = p_pi_both,
       height = 6, width = 5.5, units = "in", 
       device = "tiff", dpi = 300,
       compression = "lzw", type = "cairo")

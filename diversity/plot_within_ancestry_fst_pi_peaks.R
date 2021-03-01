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
# fst_peaks_file = "diversity/results/fst/HILO_MAIZE55/Ne10000_yesBoot/HOMOZYG/summary_pop_pairs_fst.mexicana_ancestry.peaks.allChr.txt"
fst_genomewide_file = snakemake@input[["fst_genomewide"]]
# fst_genomewide_file = "diversity/results/fst/HILO_MAIZE55/Ne10000_yesBoot/HOMOZYG/summary_pop_pairs_fst.allChr.txt"
pi_peaks_file = snakemake@input[["pi_peaks"]]
# pi_peaks_file = "diversity/results/pi/HILO_MAIZE55/Ne10000_yesBoot/HOMOZYG/summary_pop_pi.mexicana_ancestry.peaks.allChr.txt"
pi_genomewide_file = snakemake@input[["pi_genomewide"]]
# pi_genomewide_file = "diversity/results/pi/HILO_MAIZE55/Ne10000_yesBoot/HOMOZYG/summary_pop_pi.allChr.txt"
load(snakemake@input[["meta"]])
# load("samples/HILO_MAIZE55_meta.RData")
png_pi_mexicana_anc = snakemake@output[["png_pi_mex"]]
# png_pi_mexicana_anc = "diversity/plots/HILO_MAIZE55/Ne10000_yesBoot/pi_within_mexicana_ancestry_peaks.png"
png_pi_maize_anc = snakemake@output[["png_pi_maize"]]
# png_pi_maize_anc = "diversity/plots/HILO_MAIZE55/Ne10000_yesBoot/pi_within_maize_ancestry.png"
png_fst_mexicana_anc = snakemake@output[["png_fst_mex"]]
# png_fst_mexicana_anc = "diversity/plots/HILO_MAIZE55/Ne10000_yesBoot/local_fst_within_mexicana_ancestry_peaks.png"


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
                  into = c("v1", "v2", "v3", "v4", "v5", "v6", "ancestry", "pops"), 
                  sep = "/") %>%
  tidyr::separate(data = .,
                  col = pops,
                  into = c("pop1", "pop2", "v7", "v8", "v9"),
                  sep = "[.]") %>%
  dplyr::mutate(genomic_region = "genomewide") %>%
  dplyr::select(ancestry, pop1, pop2, genomic_region, fst)

fst_peaks <- read.table(fst_peaks_file, sep = " ", header = F) %>%
  data.table::setnames(c("file_name", "fst")) %>%
  tidyr::separate(data = .,
                  col = file_name, 
                  into = c("v1", "v2", "v3", "v4", "v5", "v6", "ancestry", "pops"), 
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
                  into = c("v1", "v2", "v3", "v4", "v5", "v6", "ancestry", "label"), 
                  sep = "/") %>%
  tidyr::separate(data = .,
                  col = label,
                  into = c("pop", "v7", "v8", "v9"),
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
                  into = c("v1", "v2", "v3", "v4", "v5", "v6", "ancestry", "label"), 
                  sep = "/") %>%
  tidyr::separate(data = .,
                  col = label,
                  into = c("pop", "genomic_region", "outliers_pop", "v7", "v8", "v9"),
                  sep = "[.]") %>%
  dplyr::select(ancestry, pop, genomic_region, pi)

pi = bind_rows(pi_genomewide, pi_peaks) %>%
  left_join(., meta_pops, by = "pop")


# ---------- plot data ---------- #
shapes_peaks <- c(17, 16, 0)
names(shapes_peaks) <- c("1pop", "4pop", "genomewide")
labels_peaks_maize <- c("1 population\npeaks in maize", "4+ population\npeaks in maize", "all ancestry\ngenomewide")
labels_peaks_mex <- c("1 population\npeaks in mexicana", "4+ population\npeaks in mexicana", "all ancestry\ngenomewide")


p_fst_mex <- fst %>%
  filter(LOCALITY.pop1 == LOCALITY.pop2 & 
           ancestry == "mexicana" &
           zea.pop1 == "maize") %>%
  mutate(LOCALITY.pop1 = reorder(LOCALITY.pop1, ELEVATION.pop1)) %>%
  ggplot(aes(x = LOCALITY.pop1, shape = genomic_region, y = fst)) +
  geom_point() +
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


p_pi_mex <- pi %>%
  filter(ancestry == "mexicana") %>%
  mutate(LOCALITY = reorder(LOCALITY, ELEVATION)) %>%
  ggplot(aes(x = LOCALITY, color = zea, shape = genomic_region, y = pi)) +
  geom_point() +
  coord_flip() +
  theme_light() +
  scale_shape_manual(values = shapes_peaks,
                     labels = labels_peaks_maize) +
  scale_color_manual(values = col_maize_mex_parv) +
  labs(#title = "Diversity within mexicana ancestry", 
       x = "Location", 
       y = expression(pi), 
       color = NULL, 
       shape = NULL) +
  theme(legend.key.size = unit(10, "mm"))

ggsave(filename = png_pi_mexicana_anc,
       plot = p_pi_mex,
       height = 3, width = 5.5, units = "in", 
       device = "png", dpi = 300)


p_pi_maize <- pi %>%
  filter(ancestry == "maize") %>%
  mutate(LOCALITY = reorder(LOCALITY, ELEVATION)) %>%
  ggplot(aes(x = LOCALITY, color = zea, shape = genomic_region, y = pi)) +
  geom_point() +
  coord_flip() +
  theme_light() +
  scale_shape_manual(values = shapes_peaks[3],
                     labels = labels_peaks_mex[3]) +
  scale_color_manual(values = col_maize_mex_parv) +
  labs(#title = "Diversity within maize ancestry", 
       x = "Location", 
       y = expression(pi), 
       color = NULL, 
       shape = NULL) +
  theme(legend.key.size = unit(10, "mm"))

ggsave(filename = png_pi_maize_anc,
       plot = p_pi_maize,
       height = 3, width = 5.5, units = "in", 
       device = "png", dpi = 300)

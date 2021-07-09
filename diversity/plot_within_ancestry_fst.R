#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

# script to summarise pairwise differentiation (fst) between populations
# for homozygous ancestry windows across the genome

# get snakemake variables
input_file = snakemake@input[["fst"]]
# input_file = "diversity/results/fst/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/HOMOZYG/summary_pop_pairs_fst.allChr.txt"
load(snakemake@input[["meta"]])
# load("samples/HILO_MAIZE55_PARV50_meta.RData")
png_heatmap_both = snakemake@output[["png_heatmap_both"]]
# png_heatmap_both = "diversity/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_fst_within_maize_or_mexicana_ancestry_genomewide_heatmap_both.png"
png_heatmap_both_lzw = snakemake@output[["png_heatmap_both_lzw"]]
# png_heatmap_both_lzw = "../hilo_manuscript/figures_main/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_fst_within_maize_or_mexicana_ancestry_genomewide_heatmap_both.tif"
png_heatmap_parv = snakemake@output[["png_heatmap_parv"]]
# png_heatmap_parv = "diversity/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_fst_within_parviglumis_ancestry_genomewide_heatmap.png"
png_heatmap_parv_lzw = snakemake@output[["png_heatmap_parv_lzw"]]
# png_heatmap_parv_lzw = "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_fst_within_parviglumis_ancestry_genomewide_heatmap.tif"

meta_pops = meta %>%
  dplyr::select(popN, zea, symp_allo, group, LOCALITY, ELEVATION, LAT, LON) %>%
  dplyr::distinct() %>%
  dplyr::mutate(., pop = paste0("pop", popN)) %>%
  arrange(zea, ELEVATION) %>%
  filter(symp_allo == "sympatric")

fst0 <- read.table(input_file, sep = " ", header = F) %>%
  data.table::setnames(c("file_name", "fst")) %>%
  tidyr::separate(data = .,
                  col = file_name, 
                  into = c("v1", "v2", "v3", "v4", "v5", "v6", "v7", "ancestry", "pops"), 
                  sep = "/") %>%
  tidyr::separate(data = .,
                  col = pops,
                  into = c("pop1", "pop2", "v7", "v8", "v9"),
                  sep = "[.]") %>%
  dplyr::select(ancestry, pop1, pop2, fst)

# add in metadata for both populations being compared,
# and reverse comparison so a-b and b-a are both in dataset
fst <- bind_rows(dplyr::mutate(fst0, comparison_order = 1), # because fst is symmetric, add in all 'reverse' comparisons
            dplyr::mutate(fst0, comparison_order = 2) %>%
            dplyr::rename(., pop1 = pop2, pop2 = pop1)) %>%
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

# plot the fst distance matrix as a heatmap:
fst_mexicana = filter(fst, ancestry == "mexicana")
fst_maize = filter(fst, ancestry == "maize")
fst_parv = filter(fst, ancestry == "parv")


# tiles but for parviglumis ancestry:
p_heatmap_parv <- fst_parv %>%
  arrange(., zea.pop1, ELEVATION.pop1) %>%
  mutate(INDEX.pop1 = 1:nrow(.)) %>%
  arrange(., zea.pop2, ELEVATION.pop2) %>%
  mutate(INDEX.pop2 = 1:nrow(.)) %>%
  mutate(zea_loc1 = paste(LOCALITY.pop1, zea.pop1),
         zea_loc2 = paste(LOCALITY.pop2, zea.pop2),
         zea_loc1 = reorder(zea_loc1, INDEX.pop1),
         zea_loc2 = reorder(zea_loc2, INDEX.pop2))%>%
  ggplot(data = ., 
         aes(x = zea_loc1, 
             y = zea_loc2, 
             fill = fst)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_viridis(option = "magma", direction = -1, limits = c(-0.025, 0.5)) +
  labs(x = "population 1", y = "population 2", fill = expression(F[ST])) +
  #ggtitle("fst within parviglumis ancestry") +
  coord_fixed() +
  geom_point(data = fst_parv %>%
               arrange(., zea.pop1, ELEVATION.pop1) %>%
               mutate(INDEX.pop1 = 1:nrow(.)) %>%
               arrange(., zea.pop2, ELEVATION.pop2) %>%
               mutate(INDEX.pop2 = 1:nrow(.)) %>%
               mutate(zea_loc1 = paste(LOCALITY.pop1, zea.pop1),
                      zea_loc2 = paste(LOCALITY.pop2, zea.pop2),
                      zea_loc1 = reorder(zea_loc1, INDEX.pop1),
                      zea_loc2 = reorder(zea_loc2, INDEX.pop2)) %>%
               mutate(comparison_type = ifelse(LOCALITY.pop1 == LOCALITY.pop2, "same location", "different location")) %>%
               filter(comparison_type == "same location"),
    aes(color = comparison_type),
             size = 0.5,
             pch = 19) +
  scale_color_manual(values = "white", labels = "sympatric\npopulation pair") +
  theme(legend.key = element_rect(fill = "darkgrey", color = NA),
        plot.title = element_blank()) +
  labs(color = NULL, fill = expression(F[ST]))
# p_heatmap_parv
ggsave(filename = png_heatmap_parv,
       plot = p_heatmap_parv,
       height = 4.75, width = 6, 
       units = "in", device = "png", dpi = 300)
ggsave(filename = png_heatmap_parv_lzw,
       plot = p_heatmap_parv,
       height = 4.75, width = 6, 
       units = "in", device = "tiff", dpi = 300,
       compression = "lzw", type = "cairo")

# plot both maize and mexicana ancestries with one heatmap
p_heatmap_both <- fst_maize %>%
  arrange(., zea.pop1, ELEVATION.pop1) %>%
  mutate(INDEX.pop1 = 1:nrow(.)) %>%
  arrange(., zea.pop2, ELEVATION.pop2) %>%
  mutate(INDEX.pop2 = 1:nrow(.)) %>%
  mutate(zea_loc1 = paste(LOCALITY.pop1, zea.pop1),
         zea_loc2 = paste(LOCALITY.pop2, zea.pop2),
         zea_loc1 = reorder(zea_loc1, INDEX.pop1),
         zea_loc2 = reorder(zea_loc2, INDEX.pop2)) %>%
  filter(INDEX.pop1 < INDEX.pop2 | INDEX.pop1 == INDEX.pop2) %>%
  bind_rows(.,
            fst_mexicana %>%
              arrange(., zea.pop1, ELEVATION.pop1) %>%
              mutate(INDEX.pop1 = 1:nrow(.)) %>%
              arrange(., zea.pop2, ELEVATION.pop2) %>%
              mutate(INDEX.pop2 = 1:nrow(.)) %>%
              mutate(zea_loc1 = paste(LOCALITY.pop1, zea.pop1),
                     zea_loc2 = paste(LOCALITY.pop2, zea.pop2),
                     zea_loc1 = reorder(zea_loc1, INDEX.pop1),
                     zea_loc2 = reorder(zea_loc2, INDEX.pop2)) %>%
              filter(INDEX.pop1 > INDEX.pop2 | INDEX.pop1 == INDEX.pop2)) %>%
  mutate( comparison_type = ifelse(LOCALITY.pop1 == LOCALITY.pop2, "same location", "different location")) %>%
  ggplot(., 
         aes(x = zea_loc1, 
             y = zea_loc2, 
             fill = fst)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_viridis(option = "magma", direction = -1, limits = c(-0.025, 0.5)) +
  labs(x = "population 1", y = "population 2") +
  #ggtitle("fst within maize (top left)\nor mexicana (bottom right) ancestry") +
  coord_fixed() +
  geom_point(data = fst %>%
                           arrange(., zea.pop1, ELEVATION.pop1) %>%
                           mutate(INDEX.pop1 = 1:nrow(.)) %>%
                           arrange(., zea.pop2, ELEVATION.pop2) %>%
                           mutate(INDEX.pop2 = 1:nrow(.)) %>%
                           mutate(zea_loc1 = paste(LOCALITY.pop1, zea.pop1),
                                  zea_loc2 = paste(LOCALITY.pop2, zea.pop2),
                                  zea_loc1 = reorder(zea_loc1, INDEX.pop1),
                                  zea_loc2 = reorder(zea_loc2, INDEX.pop2)) %>%
               mutate(comparison_type = ifelse(LOCALITY.pop1 == LOCALITY.pop2, "same location", "different location")) %>%
               filter(comparison_type == "same location"),
             aes(color = comparison_type),
             size = 0.5,
             pch = 19) +
  scale_color_manual(values = "white", labels = "sympatric\npopulation pair") +
  theme(legend.key = element_rect(fill = "darkgrey", color = NA),
        plot.title = element_blank()) +
  labs(color = NULL, fill = expression(F[ST]))
# note: one fst value is very small but negative
# p_heatmap_both

ggsave(filename = png_heatmap_both,
       plot = p_heatmap_both,
       height = 4.75, width = 6, 
       units = "in", device = "png", dpi = 300)

ggsave(filename = png_heatmap_both_lzw,
       plot = p_heatmap_both,
       height = 4.75, width = 6, 
       units = "in", device = "tiff", dpi = 300,
       compression = "lzw", type = "cairo")

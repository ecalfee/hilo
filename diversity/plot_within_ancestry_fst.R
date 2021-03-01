#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

# script to summarise pairwise differentiation (fst) between populations
# for homozygous ancestry windows across the genome

# get snakemake variables
input_file = snakemake@input[["fst"]]
# input_file = "diversity/results/fst/HILO_MAIZE55/Ne10000_yesBoot/HOMOZYG/summary_pop_pairs_fst.allChr.txt"
load(snakemake@input[["meta"]])
# load("samples/HILO_MAIZE55_meta.RData")
png_points_maize = snakemake@output[["png_points_maize"]]
# png_points_maize = "diversity/plots/HILO_MAIZE55/Ne10000_yesBoot/fst_within_maize_ancestry_genomewide_points.png"
png_points_mexicana = snakemake@output[["png_points_mexicana"]]
# png_points_mexicana = "diversity/plots/HILO_MAIZE55/Ne10000_yesBoot/fst_within_mexicana_ancestry_genomewide_points.png"
png_heatmap_maize = snakemake@output[["png_heatmap_maize"]]
# png_heatmap_maize = "diversity/plots/HILO_MAIZE55/Ne10000_yesBoot/fst_within_maize_ancestry_genomewide_heatmap.png"
png_heatmap_mexicana = snakemake@output[["png_heatmap_mexicana"]]
# png_heatmap_mexicana = "diversity/plots/HILO_MAIZE55/Ne10000_yesBoot/fst_within_mexicana_ancestry_genomewide_heatmap.png"
png_heatmap_both = snakemake@output[["png_heatmap_both"]]
# png_heatmap_both = "diversity/plots/HILO_MAIZE55/Ne10000_yesBoot/fst_within_maize_or_mexicana_ancestry_genomewide_heatmap_both.png"


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
                  into = c("v1", "v2", "v3", "v4", "v5", "v6", "ancestry", "pops"), 
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

# plot fst within mexicana ancestry
p_points_mexicana <- fst %>%
  dplyr::filter(ancestry == "mexicana") %>%
  mutate(symp_allo = ifelse(LOCALITY.pop1==LOCALITY.pop2, "sympatric", "allopatric")) %>%
  mutate(comparison_type = factor(paste(zea.pop1, zea.pop2, sep = "-"),
                                  ordered = T,
                                  levels = c("maize-maize", "mexicana-mexicana", "maize-mexicana", "mexicana-maize")),
         LOCALITY.pop1 = reorder(LOCALITY.pop1, ELEVATION.pop1),
         LOCALITY.pop2 = reorder(LOCALITY.pop2, -ELEVATION.pop2)) %>%
  ggplot(., aes(x = LOCALITY.pop1, 
                y = fst,
                color = LOCALITY.pop2,
                shape = symp_allo,
                size = symp_allo)) +
  geom_jitter(width = .2) +
  facet_wrap(~comparison_type, nrow = 2) +
  coord_flip() +
  theme_light() +
  labs(x = "Population 1", color = "Population 2", y = expression(F[ST]),
       shape = "Comparison Type", size = "Comparison Type") +
  scale_shape_manual(values = c(19, 1)) +
  scale_size_manual(values = c(1, 2)) +
  scale_color_viridis_d(direction = 1, option = "viridis")# +
  #ggtitle("Fst within mexicana ancestry tracts (genomewide)")

ggsave(file = png_points_mexicana,
       plot = p_points_mexicana,
       height = 7, width = 7.5, 
       units = "in", device = "png", dpi = 300)

# plot fst within maize ancestry
p_points_maize <- fst %>%
  dplyr::filter(ancestry == "maize") %>%
  mutate(symp_allo = ifelse(LOCALITY.pop1==LOCALITY.pop2, "sympatric", "allopatric")) %>%
  mutate(comparison_type = factor(paste(zea.pop1, zea.pop2, sep = "-"),
                                  ordered = T,
                                  levels = c("maize-maize", "mexicana-mexicana", "maize-mexicana", "mexicana-maize")),
         LOCALITY.pop1 = reorder(LOCALITY.pop1, ELEVATION.pop1),
         LOCALITY.pop2 = reorder(LOCALITY.pop2, -ELEVATION.pop2)) %>%
  ggplot(., aes(x = LOCALITY.pop1, 
                y = fst,
                color = LOCALITY.pop2,
                shape = symp_allo,
                size = symp_allo)) +
  geom_jitter(width = .2) +
  facet_wrap(~comparison_type, nrow = 2) +
  coord_flip() +
  theme_light() +
  labs(x = "Population 1", color = "Population 2", y = expression(F[ST]),
       shape = "Comparison Type", size = "Comparison Type") +
  scale_shape_manual(values = c(19, 1)) +
  scale_size_manual(values = c(1, 2)) +
  scale_color_viridis_d(direction = 1, option = "viridis")# +
  #ggtitle("Fst within maize ancestry tracts (genomewide)")

ggsave(file = png_points_maize,
       plot = p_points_maize,
       height = 7, width = 7.5, 
       units = "in", device = "png", dpi = 300)


# plot the fst distance matrix as a heatmap:
fst_mexicana = filter(fst, ancestry == "mexicana")
fst_maize = filter(fst, ancestry == "maize")

p_heatmap_mexicana <- fst_mexicana %>%
  arrange(., zea.pop1, ELEVATION.pop1) %>%
  mutate(INDEX.pop1 = 1:nrow(.)) %>%
  arrange(., zea.pop2, ELEVATION.pop2) %>%
  mutate(INDEX.pop2 = 1:nrow(.)) %>%
  mutate(zea_loc1 = paste(LOCALITY.pop1, zea.pop1),
         zea_loc2 = paste(LOCALITY.pop2, zea.pop2),
         zea_loc1 = reorder(zea_loc1, INDEX.pop1),
         zea_loc2 = reorder(zea_loc2, INDEX.pop2)) %>%
  ggplot(., 
         aes(x = zea_loc1, 
             y = zea_loc2, 
             fill = fst)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_viridis(option = "magma", direction = -1, limits = c(-0.01, 0.5)) +
  labs(x = "population 1", y = "population 2", fill = expression(F[ST])) +
  #ggtitle("fst within mexicana ancestry") +
  coord_fixed()

ggsave(file = png_heatmap_mexicana,
       plot = p_heatmap_mexicana,
       height = 6, width = 6, 
       units = "in", device = "png", dpi = 300)

# tiles but for maize ancestry:
p_heatmap_maize <- fst_maize %>%
  arrange(., zea.pop1, ELEVATION.pop1) %>%
  mutate(INDEX.pop1 = 1:nrow(.)) %>%
  arrange(., zea.pop2, ELEVATION.pop2) %>%
  mutate(INDEX.pop2 = 1:nrow(.)) %>%
  mutate(zea_loc1 = paste(LOCALITY.pop1, zea.pop1),
         zea_loc2 = paste(LOCALITY.pop2, zea.pop2),
         zea_loc1 = reorder(zea_loc1, INDEX.pop1),
         zea_loc2 = reorder(zea_loc2, INDEX.pop2)) %>%
  ggplot(., 
         aes(x = zea_loc1, 
             y = zea_loc2, 
             fill = fst)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_viridis(option = "magma", direction = -1, limits = c(-0.01, 0.5)) +
  labs(x = "population 1", y = "population 2", fill = expression(F[ST])) +
  #ggtitle("fst within maize ancestry") +
  coord_fixed()
# note: one fst value is very small but negative

ggsave(file = png_heatmap_maize,
       plot = p_heatmap_maize,
       height = 6, width = 6, 
       units = "in", device = "png", dpi = 300)


# plot both ancestries with one heatmap
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
  scale_fill_viridis(option = "magma", direction = -1, limits = c(-0.01, 0.5)) +
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
#p_heatmap_both

ggsave(file = png_heatmap_both,
       plot = p_heatmap_both,
       height = 4.75, width = 6, 
       units = "in", device = "png", dpi = 300)

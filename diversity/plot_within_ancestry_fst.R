#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(ape)
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
png_tree_maize = snakemake@output[["png_tree_maize"]]
# png_tree_maize = "diversity/plots/HILO_MAIZE55/Ne10000_yesBoot/fst_within_maize_ancestry_genomewide_tree.png"
png_tree_mexicana = snakemake@output[["png_tree_mexicana"]]
# png_tree_mexicana = "diversity/plots/HILO_MAIZE55/Ne10000_yesBoot/fst_within_mexicana_ancestry_genomewide_tree.png"
png_tree_fit_maize = snakemake@output[["png_tree_fit_maize"]]
# png_tree_fit_maize = "diversity/plots/HILO_MAIZE55/Ne10000_yesBoot/fst_within_maize_ancestry_genomewide_tree_fit.png"
png_tree_fit_mexicana = snakemake@output[["png_tree_fit_mexicana"]]
# png_tree_fit_mexicana = "diversity/plots/HILO_MAIZE55/Ne10000_yesBoot/fst_within_mexicana_ancestry_genomewide_tree_fit.png"
png_heatmap_maize = snakemake@output[["png_heatmap_maize"]]
# png_heatmap_maize = "diversity/plots/HILO_MAIZE55/Ne10000_yesBoot/fst_within_maize_ancestry_genomewide_heatmap.png"
png_heatmap_mexicana = snakemake@output[["png_heatmap_mexicana"]]
# png_heatmap_mexicana = "diversity/plots/HILO_MAIZE55/Ne10000_yesBoot/fst_within_mexicana_ancestry_genomewide_heatmap.png"

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
  labs(x = "Population 1", color = "Population 2", 
       shape = "Comparison Type", size = "Comparison Type") +
  scale_shape_manual(values = c(19, 1)) +
  scale_size_manual(values = c(1, 2)) +
  scale_color_viridis_d(direction = 1, option = "viridis") +
  ggtitle("Fst within mexicana ancestry tracts (genomewide)")

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
  labs(x = "Population 1", color = "Population 2", 
       shape = "Comparison Type", size = "Comparison Type") +
  scale_shape_manual(values = c(19, 1)) +
  scale_size_manual(values = c(1, 2)) +
  scale_color_viridis_d(direction = 1, option = "viridis") +
  ggtitle("Fst within maize ancestry tracts (genomewide)")

ggsave(file = png_points_maize,
       plot = p_points_maize,
       height = 7, width = 7.5, 
       units = "in", device = "png", dpi = 300)

# make a neighbor-joining tree for fst

# first, make a distance matrix using pairwise fst
dist_matrix_mexicana <- matrix(0, 
                      nrow(meta_pops),
                      nrow(meta_pops),
                      dimnames = list(meta_pops$pop, meta_pops$pop))
fst_mexicana = filter(fst, ancestry == "mexicana")
for (r in 1:nrow(fst_mexicana)){
  dist_matrix_mexicana[fst_mexicana$pop1[r], fst_mexicana$pop2[r]] <- fst$fst[r]
}
tree_mexicana <- nj(X = dist_matrix_mexicana)

# colors for tree tips
col_locality <- viridis_pal(direction = -1, option = "viridis")(14)
names(col_locality) <- unique(meta_pops$LOCALITY)

# plot the tree mexicana ancestry
png(filename = png_tree_mexicana,
    width = 6, height = 4, 
    units = "in", res = 300)
plot(tree_mexicana, show.tip = F, 
     cex = 0.6, x.lim = c(0, 0.4))
title("Mexicana ancestry - NJ Tree (Fst)")
tiplabels(text = paste(meta_pops$LOCALITY, meta_pops$zea),#tree_mexicana$tip.label, 
          col = col_locality[meta_pops$LOCALITY],
          bg = NULL,
          frame = "none",
          adj = -0.25,
          cex = 0.6)
ape::axisPhylo()
dev.off()

# how good is nj tree fit? ok
x_mexicana <- as.vector(dist_matrix_mexicana)
y_mexicana <- as.vector(cophenetic(tree_mexicana))
png(filename = png_tree_fit_mexicana,
    width = 5, height = 5, 
    units = "in", res = 300)
plot(x_mexicana, y_mexicana, 
     xlim = c(0, max(x_mexicana, y_mexicana)),
     ylim = c(0, max(x_mexicana, y_mexicana)),
     xlab = "original pairwise fst", 
     ylab = "pairwise distances on the NJ tree",
     main = paste("Fit of NJ tree to fst within mexicana ancestry\nrho =", 
                 round(cor(x_mexicana, y_mexicana, 
                           method = "pearson"), 2)))
abline(a = 0, b = 1, col = "blue")
#abline(lm(y_mexicana ~ x_mexicana), col = "orange")
dev.off()



# now make the tree for maize ancestry:
# first, make a distance matrix using pairwise fst
dist_matrix_maize <- matrix(0, 
                               nrow(meta_pops),
                               nrow(meta_pops),
                               dimnames = list(meta_pops$pop, meta_pops$pop))
fst_maize = filter(fst, ancestry == "maize")
for (r in 1:nrow(fst_maize)){
  dist_matrix_maize[fst_maize$pop1[r], fst_maize$pop2[r]] <- fst$fst[r]
}
tree_maize <- nj(X = dist_matrix_maize)

# plot the tree
png(filename = png_tree_maize,
    width = 6, height = 4, 
    units = "in", res = 300)
plot(tree_maize, show.tip = F, 
     cex = 0.6, x.lim = c(0, 0.4))
title("Maize ancestry - NJ Tree (Fst)")
tiplabels(text = paste(meta_pops$LOCALITY, meta_pops$zea),#tree_maize$tip.label, 
          col = col_locality[meta_pops$LOCALITY],
          bg = NULL,
          frame = "none",
          adj = -0.25,
          cex = 0.6)
ape::axisPhylo()
dev.off()

# how good is nj tree fit? ok
x_maize <- as.vector(dist_matrix_maize)
y_maize <- as.vector(cophenetic(tree_maize))
png(filename = png_tree_fit_maize,
    width = 5, height = 5, 
    units = "in", res = 300)
plot(x_maize, y_maize, 
     xlim = c(0, max(x_maize, y_maize)),
     ylim = c(0, max(x_maize, y_maize)),
     xlab = "original pairwise fst", 
     ylab = "pairwise distances on the NJ tree",
     main = paste("Fit of NJ tree to fst within maize ancestry\nrho =", 
                  round(cor(x_maize, y_maize, 
                            method = "pearson"), 2)))
abline(a = 0, b = 1, col = "blue")
#abline(lm(y_maize ~ x_maize), col = "orange")
dev.off()


# rather than a tree, draw the distance matrix:
p_heatmap_mexicana <- fst_mexicana %>%
  arrange(., zea.pop1, ELEVATION.pop1) %>%
  mutate(INDEX.pop1 = 1:nrow(.)) %>%
  arrange(., zea.pop2, ELEVATION.pop2) %>%
  mutate(INDEX.pop2 = 1:nrow(.)) %>%
  mutate(zea_loc1 = paste(zea.pop1, LOCALITY.pop1),
         zea_loc2 = paste(zea.pop2, LOCALITY.pop2),
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
  labs(x = "population 1", y = "population 2") +
  ggtitle("fst within mexicana ancestry") +
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
  mutate(zea_loc1 = paste(zea.pop1, LOCALITY.pop1),
         zea_loc2 = paste(zea.pop2, LOCALITY.pop2),
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
  labs(x = "population 1", y = "population 2") +
  ggtitle("fst within maize ancestry") +
  coord_fixed()
# note: one fst value is very small but negative

ggsave(file = png_heatmap_maize,
       plot = p_heatmap_maize,
       height = 6, width = 6, 
       units = "in", device = "png", dpi = 300)

# test if 'local sympatric' pops have lower fst? maybe use ~local adaptation test.
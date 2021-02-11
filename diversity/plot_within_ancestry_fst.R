#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)

# script to summarise diversity (pi) for homozygous ancestry windows within high introgression peaks vs. not

# get snakemake variables
input_file = snakemake@input[["fst"]]
# input_file = "diversity/results/fst/HILO_MAIZE55/Ne10000_yesBoot/HOMOZYG/summary_pop_pairs_fst.allChr.txt"
load(snakemake@input[["meta"]])
# load("samples/HILO_MAIZE55_meta.RData")
png_out = snakemake@output[["png_out"]]
# png_out = "diversity/plots/HILO_MAIZE55/Ne10000_yesBoot/fst_within_ancestry_genomewide.png"

meta_pops = meta %>%
  dplyr::select(popN, zea, symp_allo, group, LOCALITY, ELEVATION, LAT, LON) %>%
  dplyr::distinct() %>%
  dplyr::mutate(., pop = paste0("pop", popN))

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
fst %>%
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


# plot fst within maize ancestry
fst %>%
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

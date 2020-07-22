#!/usr/bin/env Rscript
# working directory is hilo/
# calculation/plotting global ancestry for low-coverage maize/mexicana individuals 
# in HILO adaptation project
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(viridis)

# load variables from snakemake
# get output plot filenames
png_elev = snakemake@output[["png_elev"]]
# png_elev = "global_ancestry/plots/lm_mexicana_by_pop_elevation_K2.png"
png_structure = snakemake@output[["png_structure"]]
# png_structure = "global_ancestry/plots/structure_K2.png"

# get colors for plot
# source("colors.R")
source(snakemake@input[["colors"]])

# get NGSAdmix results for K = 2, d_admix2
load(snakemake@input[["k2"]])
# load("global_ancestry/results/NGSAdmix/HILO_MAIZE55/K2_alphas_by_ind.RData")



# make STRUCTURE-like ancestry plots:
p_structure <- d_admix2 %>%
  mutate(ELEVATION = ifelse(group == "allopatric_maize", 0, ELEVATION)) %>% # add dummy elevation for allopatric maize (from Palmar Chico)
  arrange(., ELEVATION, popN) %>%
  mutate(sample = 1:nrow(.)) %>%
  tidyr::gather(., "ancestry", "p", c("mexicana", "maize")) %>%
  ggplot(., aes(fill = ancestry, y = p, x = sample)) +
  geom_bar(stat = "identity", position = "fill") + 
  facet_grid(symp_allo ~ zea) +
  #facet_grid(zea ~ LOCALITY) +
  scale_fill_manual(values = col_maize_mex_parv) +
  ggtitle("K=2 NGSAdmix results") +
  labs(fill = "Ancestry", x = "Sample") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
#plot(p_structure)
ggsave(png_structure, 
        plot = p_structure, 
        device = "png", 
        width = 7.5, height = 4, units = "in",
        dpi = 300)

lm_maize = filter(d_admix2, group == "sympatric_maize") %>%
  with(., lm(mexicana ~ ELEVATION))
lm_mex = filter(d_admix2, group == "sympatric_mexicana") %>%
  with(., lm(mexicana ~ ELEVATION))
#glance(lm_maize)
#summary(lm_mex)
#glance(lm_mex)
#summary(lm_maize)

# print linear model output to a file?

p_symp_elev <- d_admix2 %>%
  filter(., symp_allo == "sympatric") %>%
  arrange(., ELEVATION) %>%
  mutate(LOCALITY = factor(LOCALITY, ordered = T, levels = unique(.$LOCALITY))) %>%
  ggplot(., aes(x = ELEVATION, 
                y = mexicana, 
                color = LOCALITY,
                shape = zea)) +
  geom_point(alpha = 0.75, size = 2) +
  ylab("Proportion mexicana ancestry") +
  xlab("Elevation (m)") +
  geom_abline(intercept = c(coef(lm_mex)[1], coef(lm_maize)[1]),
              slope = c(coef(lm_mex)[2], coef(lm_maize)[2])) +
  ggtitle("Clines in mexicana ancestry across elevation") +
  labs(color = "Location", shape = "Subspecies") +
  theme_classic() +
  scale_color_viridis_d(direction = -1, option = "viridis")
#p_symp_elev
ggsave(filename = png_elev,
       plot = p_symp_elev,
       device = "png", height = 6, width = 7.5, units = "in", dpi = 300)

#!/usr/bin/env Rscript
# working directory is hilo/
# calculation/plotting global ancestry for low-coverage maize/mexicana individuals 
# in HILO adaptation project
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(viridis)
library(xtable)
library(grid)
library(gridExtra)

# load variables from snakemake
# get output plot and table filenames
png_elev = snakemake@output[["png_elev"]]
# png_elev = "global_ancestry/plots/lm_mexicana_by_pop_elevation_K2.png"
png_elev_symp_allo = snakemake@output[["png_elev_symp_allo"]]
# png_elev_symp_allo = "global_ancestry/plots/lm_mexicana_by_pop_elevation_K2_symp_allo.png"
png_structure = snakemake@output[["png_structure"]]
# png_structure = "global_ancestry/plots/structure_K2.png"
lm_tex = snakemake@output[["lm_tex"]]
# lm_tex = "global_ancestry/tables/lm_elevation.tex"
png_global_anc_multi = snakemake@output[["png_global_anc_multi"]]
# png_global_anc_multi = "global_ancestry/plots/global_anc_multi.png"
png_global_anc_multi_lwz = snakemake@output[["png_global_anc_multi_lzw"]]
# png_global_anc_multi_lzw = "../hilo_manuscript/figures_main/global_anc_multi.tif"

# some plot sizing constants
smallPointSize = 0.25
smallLabel = 7
smallLegendText = 6
smallLegendSize = 2

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
        width = 5.4, height = 3, units = "in",
        dpi = 300)

lm_maize = filter(d_admix2, group == "sympatric_maize") %>%
  mutate(elevation_km = ELEVATION/1000) %>% # convert meters to km
  with(., lm(mexicana ~ elevation_km))
lm_mex = filter(d_admix2, group == "sympatric_mexicana") %>%
  mutate(elevation_km = ELEVATION/1000) %>% # convert meters to km
  with(., lm(mexicana ~ elevation_km))
#glance(lm_maize)
#summary(lm_mex)
#glance(lm_mex)
#summary(lm_maize)

lm_table <- bind_rows(mutate(broom::tidy(lm_maize), subspecies = "maize"),
                      mutate(broom::tidy(lm_mex), subspecies = "mexicana")) %>%
  dplyr::select(subspecies, term, estimate, std.error, statistic, p.value) %>%
  mutate(term = ifelse(term == "(Intercept)", "intercept", ifelse(term == "elevation_km", "elevation (km)", term)))

# print linear model output to a file
print(xtable(lm_table,
             type = "latex",
             digits = c(1, 1, 1, 4, 4, 1, -2),
             latex.environments = NULL),
      include.rownames = F,
      file = lm_tex)


p_symp_elev <- d_admix2 %>%
  filter(., symp_allo == "sympatric") %>%
  arrange(., ELEVATION) %>%
  mutate(zea = reorder(zea, desc(ELEVATION))) %>%
  mutate(LOCALITY = factor(LOCALITY, ordered = T, levels = unique(.$LOCALITY))) %>%
  ggplot(., aes(x = ELEVATION, 
                y = mexicana, 
                color = LOCALITY,
                shape = zea)) +
  geom_point(alpha = 0.75, size = 2) +
  ylab("Proportion mexicana ancestry") +
  xlab("Elevation (m)") +
  geom_abline(intercept = c(coef(lm_mex)[1], coef(lm_maize)[1]),
              slope = c(coef(lm_mex)[2], coef(lm_maize)[2])/1000) + # divided by 1000 to put on meters, not km, x-axis scale
  ggtitle("Clines in mexicana ancestry across elevation") +
  labs(color = "Location", shape = "Subspecies") +
  theme_classic() +
  scale_color_viridis_d(direction = -1, option = "viridis")
# p_symp_elev
ggsave(filename = png_elev,
       plot = p_symp_elev,
       device = "png", height = 6, width = 7.5, units = "in", dpi = 300)

# allopatric mexicana/maize structure only
p_structure_allo <- d_admix2 %>%
  filter(symp_allo == "allopatric") %>%
  mutate(ELEVATION = ifelse(group == "allopatric_maize", 983, ELEVATION)) %>%
  arrange(., ELEVATION, popN) %>%
  mutate(sample = 1:nrow(.)) %>%
  mutate(sample = ifelse(LOCALITY == "Puerta Encantada", sample + 1, 
                         ifelse(LOCALITY == "Malinalco", sample + 2, 
                                ifelse(LOCALITY == "Amecameca", sample + 3, 
                                       sample)))) %>%
  tidyr::gather(., "ancestry", "p", c("mexicana", "maize")) %>%
  ggplot(., aes(fill = ancestry, y = p, x = sample)) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(values = col_maize_mex_parv, labels = c("maize ancestry", "mexicana ancestry")) +
  labs(fill = NULL, 
       y = "Proportion") +
  theme_classic() +
  theme(axis.ticks.x = element_blank()) +
  scale_x_discrete(name = element_blank(), 
                   limits = c(27.5, 5.5 + 56, 20 + 56, 37.5 + 56),
                   labels = c("Palmar Chico\n(983m)",
                              "Puerta Encantada\n(1,658m)", 
                              "Malinalco\n(1,887m)", 
                              "Amecameca\n(2,467m)")) +
  theme(plot.margin = margin(c(0,0,0,5.5)))
# p_structure_allo

# multipanel plot of global ancestry data:
p_combined <- grid.arrange(grobs = list(ggplotGrob(p_symp_elev +
                                                     guides(color = F) +
                                                     labs(shape = "Sympatric") +
                                                     scale_shape_discrete(labels = c("Sympatric\nmexicana", "Sympatric\nmaize")) +
                                                     theme(plot.title = element_blank(),
                                                           legend.position = "right",
                                                           legend.text = element_text(size = smallLabel + 2),
                                                           legend.title = element_blank(),
                                                           legend.key.size = unit(10, "mm"),
                                                           legend.margin = margin(c(0,0,0,0)))),
                                        ggplotGrob(p_structure_allo + 
                                                     labs(subtitle = " Allopatric maize                                                                                     Allopatric mexicana") +
                                                     theme(axis.text.x = element_text(size = smallLabel + 1),
                                                           legend.key.size = unit(4, "mm"),
                                                           plot.subtitle = element_text(size = smallLabel + 2),
                                                           axis.title = element_text(size = smallLabel + 2),
                                                           legend.margin = margin(c(0,0,0,0)),
                                                           legend.position = "bottom")
                                        ),
                                        textGrob(label = "A", 
                                                 x = unit(0.5, "lines"), 
                                                 y = unit(0, "lines")),
                                        textGrob(label = "B", 
                                                 x = unit(0.5, "lines"), 
                                                 y = unit(0, "lines"))
                                       ),
                                        layout_matrix = rbind(
                                          c(4),
                                          c(2),
                                          c(5),
                                          c(1)),
                           heights = c(.1, 2, .1, 6),
                           widths = c(1))

# p_combined
ggsave(png_global_anc_multi, 
       plot = p_combined, 
       device = "png", 
       width = 7.5, 
       height = 7.5, 
       units = "in",
       dpi = 300)

ggsave(filename = png_global_anc_multi_lzw,
       plot = p_combined, 
       device = "tiff", 
       width = 7.5, 
       height = 7.5, 
       units = "in",
       dpi = 300, 
       compression = "lzw", 
       type = "cairo")
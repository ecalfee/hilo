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
png_structure = snakemake@output[["png_structure"]]
# png_structure = "global_ancestry/plots/HILO_MAIZE55_PARV50_structure_K3.png"
png_structure_lzw = snakemake@output[["png_structure_lzw"]]
# png_structure_lzw = "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_structure_K3.tif"
png_elev = snakemake@output[["png_elev"]]
# png_elev = "global_ancestry/plots/HILO_MAIZE55_PARV50_lm_mexicana_by_pop_elevation_K3.png"
lm_tex = snakemake@output[["lm_tex"]]
# lm_tex = "global_ancestry/tables/HILO_MAIZE55_PARV50_lm_elevation_K3.tex"
png_global_anc_multi = snakemake@output[["png_global_anc_multi"]]
# png_global_anc_multi = "global_ancestry/plots/HILO_MAIZE55_PARV50_global_anc_multi_K3.png"
png_global_anc_multi_lzw = snakemake@output[["png_global_anc_multi_lzw"]]
# png_global_anc_multi_lzw = "../hilo_manuscript/figures_main/HILO_MAIZE55_PARV50_global_anc_multi_K3.tif"
png_elev_parv = snakemake@output[["png_elev_parv"]]
# png_elev_parv = "global_ancestry/plots/HILO_MAIZE55_PARV50_lm_parviglumis_by_pop_elevation_K3.png"
png_elev_maize = snakemake@output[["png_elev_maize"]]
# png_elev_maize = "global_ancestry/plots/HILO_MAIZE55_PARV50_lm_maize_by_pop_elevation_K3.png"
png_elev_parv_lzw = snakemake@output[["png_elev_parv_lzw"]]
# png_elev_parv_lzw = "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_lm_parviglumis_by_pop_elevation_K3.tif"
png_elev_maize_lzw = snakemake@output[["png_elev_maize_lzw"]]
# png_elev_maize_lzw = "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_lm_maize_by_pop_elevation_K3.tif"


# get colors for plot
# source("colors.R")
source(snakemake@input[["colors"]])

# get NGSAdmix results (d_admix2)
load(snakemake@input[["k3"]])
# load("global_ancestry/results/NGSAdmix/HILO_MAIZE55_PARV50/K3_alphas_by_ind.RData")

# some plot sizing constants
smallPointSize = 0.25
smallLabel = 7
smallLegendText = 6
smallLegendSize = 2

# make STRUCTURE-like ancestry plots:
p_structure_symp <- d_admix2 %>%
  filter(symp_allo == "sympatric") %>%
  arrange(., ELEVATION, popN) %>%
  group_by(zea, LOCALITY) %>%
  mutate(sample = row_number()) %>%
  ungroup() %>%
  tidyr::gather(., "ancestry", "p", c("mexicana", "maize", "parviglumis")) %>%
  mutate(ancestry = factor(ancestry, ordered = T, levels = c("maize", "parviglumis", "mexicana"))) %>%
  ggplot(., aes(fill = ancestry, y = p, x = sample)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(zea ~ reorder(paste0(LOCALITY, "\n", "(", ELEVATION, "m", ")"),
                                  ELEVATION),
             scales = "free_x", space = "free_x") +
  scale_fill_manual(values = col_maize_mex_parv) +
  labs(fill = "Ancestry", x = "Sample") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(angle = 90)) +
  ylab("Admixture proportion")

# p_structure_symp

# allopatric mexicana/maize/parviglumis structure only
p_structure_allo <- d_admix2 %>%
  filter(group %in% c("allopatric_maize", "allopatric_mexicana", "parviglumis")) %>%
  mutate(ELEVATION = ifelse(group == "allopatric_maize", 983, ELEVATION)) %>%
  mutate(ELEVATION = ifelse(group == "parviglumis", 1008, ELEVATION)) %>%
  arrange(., ELEVATION, popN) %>%
  mutate(sample = 1:nrow(.)) %>%
  mutate(sample = ifelse(group == "parviglumis", sample + 1*2,
                         ifelse(LOCALITY == "Puerta Encantada", sample + 2*2,
                                ifelse(LOCALITY == "Malinalco", sample + 3*2,
                                       ifelse(LOCALITY == "Amecameca", sample + 4*2,
                                              sample))))) %>%
  tidyr::gather(., "ancestry", "p", c("mexicana", "maize", "parviglumis")) %>%
  mutate(ancestry = factor(ancestry, ordered = T, levels = c("maize", "parviglumis", "mexicana"))) %>%
  ggplot(., aes(fill = ancestry, y = p, x = sample)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = col_maize_mex_parv[c(1,3,2)], 
                    labels = c("maize ancestry", "mexicana ancestry", "parviglumis ancestry")[c(1,3,2)]) +
  labs(fill = NULL,
       y = "Admixture Proportion") +
  theme_classic() +
  theme(axis.ticks.x = element_blank()) +
  scale_x_discrete(name = element_blank(),
                   limits = c(28, 25.5 + 57, 5.5 + 57 + 52, 21 + 57 + 52, 40 + 57 + 52),
                   #limits = c(27.5, 25 + 56, 5.5 + 56 + 51, 20 + 56 + 51, 37.5 + 56 + 51),
                   labels = c("Maize\nPalmar Chico\n(983m)",
                              "Parviglumis\nPalmar Chico\n(1008m)",
                              "Mexicana\nP. Encantada\n(1658m)",
                              "Mexicana\nMalinalco\n(1887m)",
                              "Mexicana\nAmecameca\n(2467m)")) +
  theme(plot.margin = margin(c(0,0,0,5.5)))
# p_structure_allo


# multipanel plot of global ancestry data:
p_structure <- grid.arrange(grobs = list(ggplotGrob(p_structure_allo +
                                                     theme(
                                                       axis.text.x = element_text(size = 8),
                                                       legend.position = "None")),
                                        ggplotGrob(p_structure_symp +
                                                     theme(legend.position = "None")),
                                        cowplot::get_legend(p_structure_allo +
                                                              theme(legend.position = "bottom")),
                                        textGrob(label = "A",
                                                 x = unit(0.5, "lines"),
                                                 y = unit(0, "lines")),
                                        textGrob(label = "B",
                                                 x = unit(0.5, "lines"),
                                                 y = unit(0, "lines"))
                                       ),
                                        layout_matrix = rbind(
                                          c(4),
                                          c(3),
                                          c(1),
                                          c(5),
                                          c(2)),
                           heights = c(.1, .6, 3, .1, 7),
                           widths = c(1))

# p_structure
ggsave(png_structure,
       plot = p_structure,
       device = "png",
       width = 7.5,
       height = 5.5,
       units = "in",
       dpi = 300)

# save as compressed tiff
ggsave(filename = png_structure_lzw,
       plot = p_structure,
       device = "tiff",
       width = 7.5,
       height = 5.5,
       units = "in",
       dpi = 300,
       compression = "lzw",
       type = "cairo")

# how does ancestry change across elevation? linear model ancestry ~ elevation
# mexicana ancestry ~ elevation
lm_maize = filter(d_admix2, group == "sympatric_maize") %>%
  mutate(elevation_km = ELEVATION/1000) %>% # convert meters to km
  with(., lm(mexicana ~ elevation_km))
lm_mex = filter(d_admix2, group == "sympatric_mexicana") %>%
  mutate(elevation_km = ELEVATION/1000) %>% # convert meters to km
  with(., lm(mexicana ~ elevation_km))
print("linear model: mex ~ elev (km) in sympatric maize")
glance(lm_maize)
summary(lm_maize)
print("linear model: mex ~ elev (km) in sympatric mexicana")
glance(lm_mex)
summary(lm_mex)

# parviglumis ancestry ~ elevation
lm_maize_parv = filter(d_admix2, group == "sympatric_maize") %>%
  mutate(elevation_km = ELEVATION/1000) %>% # convert meters to km
  with(., lm(parviglumis ~ elevation_km))
lm_mex_parv = filter(d_admix2, group == "sympatric_mexicana") %>%
  mutate(elevation_km = ELEVATION/1000) %>% # convert meters to km
  with(., lm(parviglumis ~ elevation_km))
print("linear model: parv ~ elev (km) in sympatric maize")
glance(lm_maize_parv)
summary(lm_maize_parv)
print("linear model: parv ~ elev (km) in sympatric mexicana")
glance(lm_mex_parv)
summary(lm_mex_parv)

# maize ancestry ~ elevation
lm_maize_maize = filter(d_admix2, group == "sympatric_maize") %>%
  mutate(elevation_km = ELEVATION/1000) %>% # convert meters to km
  with(., lm(maize ~ elevation_km))
lm_mex_maize = filter(d_admix2, group == "sympatric_mexicana") %>%
  mutate(elevation_km = ELEVATION/1000) %>% # convert meters to km
  with(., lm(maize ~ elevation_km))
print("linear model: maize ~ elev (km) in sympatric maize")
glance(lm_maize_maize)
summary(lm_maize_maize)
print("linear model: maize ~ elev (km) in sympatric mexicana")
glance(lm_mex_maize)
summary(lm_mex_maize)

# all model results
lm_table <- bind_rows(mutate(broom::tidy(lm_maize), subspecies = "maize", model = "mex ~ elev"),
                      mutate(broom::tidy(lm_mex), subspecies = "mexicana", model = "mex ~ elev"),
                      mutate(broom::tidy(lm_maize_parv), subspecies = "maize", model = "parv ~ elev"),
                      mutate(broom::tidy(lm_mex_parv), subspecies = "mexicana", model = "parv ~ elev"),
                      mutate(broom::tidy(lm_maize_maize), subspecies = "maize", model = "maize ~ elev"),
                      mutate(broom::tidy(lm_mex_maize), subspecies = "mexicana", model = "maize ~ elev")) %>%
  dplyr::select(model, subspecies, term, estimate, std.error, statistic, p.value) %>%
  mutate(term = ifelse(term == "(Intercept)", "intercept", ifelse(term == "elevation_km", "elevation (km)", term)))


# print linear model output to a file
print(xtable(lm_table,
             type = "latex",
             digits = c(1, 1, 1, 1, 4, 4, 1, -2),
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
  ylab("Mexicana ancestry proportion") +
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


# combined main figure
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
                                                     ylab("Proportion") +
                                                     theme(axis.text.x = element_text(size = smallLabel + 1),
                                                           legend.key.size = unit(4, "mm"),
                                                           plot.subtitle = element_text(size = smallLabel + 2),
                                                           axis.title = element_text(size = smallLabel + 2),
                                                           legend.margin = margin(c(0,0,0,0)),
                                                           legend.position = "top"
                                                           )
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


# plot parviglumis ancestry ~ elevation in sympatric populations
p_symp_elev_parv <- d_admix2 %>%
  filter(., symp_allo == "sympatric") %>%
  arrange(., ELEVATION) %>%
  mutate(zea = reorder(zea, desc(ELEVATION))) %>%
  mutate(LOCALITY = factor(LOCALITY, ordered = T, levels = unique(.$LOCALITY))) %>%
  ggplot(., aes(x = ELEVATION, 
                y = parviglumis, 
                color = LOCALITY,
                shape = zea)) +
  geom_point(alpha = 0.75, size = 2) +
  ylab("Parviglumis ancestry proportion") +
  xlab("Elevation (m)") +
  geom_abline(data = data.frame(
    intercept = c(coef(lm_mex_parv)[1], coef(lm_maize_parv)[1]),
    slope = c(coef(lm_mex_parv)[2], coef(lm_maize_parv)[2])/1000,
    group = c("sympatric_mexicana", "sympatric_maize"),
    stringsAsFactors = F),
    aes(intercept = intercept, slope = slope)) + # divided by 1000 to put on meters, not km, x-axis scale
  #ggtitle("Clines in parviglumis ancestry across elevation") +
  labs(color = "Location", shape = "Subspecies") +
  theme_classic() +
  scale_color_viridis_d(direction = -1, option = "viridis") +
  facet_wrap(~ group, labeller = labeller(group = zea_group_labels))
# p_symp_elev_parv
ggsave(filename = png_elev_parv,
       plot = p_symp_elev_parv,
       device = "png", height = 6, width = 7.5, 
       units = "in", dpi = 300)
ggsave(filename = png_elev_parv_lzw,
       plot = p_symp_elev_parv, 
       device = "tiff", 
       width = 6, 
       height = 7.5, 
       units = "in",
       dpi = 300, 
       compression = "lzw", 
       type = "cairo")

# plot maize ancestry ~ elevation in sympatric populations
p_symp_elev_maize <- d_admix2 %>%
  filter(., symp_allo == "sympatric") %>%
  arrange(., ELEVATION) %>%
  mutate(zea = reorder(zea, desc(ELEVATION))) %>%
  mutate(LOCALITY = factor(LOCALITY, ordered = T, levels = unique(.$LOCALITY))) %>%
  ggplot(., aes(x = ELEVATION, 
                y = maize, 
                color = LOCALITY,
                shape = zea)) +
  geom_point(alpha = 0.75, size = 2) +
  ylab("Maize ancestry proportion") +
  xlab("Elevation (m)") +
  geom_abline(data = data.frame(
    intercept = c(coef(lm_mex_maize)[1], coef(lm_maize_maize)[1]),
    slope = c(coef(lm_mex_maize)[2], coef(lm_maize_maize)[2])/1000,
    group = c("sympatric_mexicana", "sympatric_maize"),
    stringsAsFactors = F),
    aes(intercept = intercept, slope = slope)) + # divided by 1000 to put on meters, not km, x-axis scale
  #ggtitle("Clines in maize ancestry across elevation") +
  labs(color = "Location", shape = "Subspecies") +
  theme_classic() +
  scale_color_viridis_d(direction = -1, option = "viridis") +
  facet_wrap(~ group, labeller = labeller(group = zea_group_labels))
# p_symp_elev_maize
ggsave(filename = png_elev_maize,
       plot = p_symp_elev_maize,
       device = "png", height = 6, width = 7.5, 
       units = "in", dpi = 300)
ggsave(filename = png_elev_maize_lzw,
       plot = p_symp_maize_maize, 
       device = "tiff", 
       width = 6, 
       height = 7.5, 
       units = "in",
       dpi = 300, 
       compression = "lzw", 
       type = "cairo")
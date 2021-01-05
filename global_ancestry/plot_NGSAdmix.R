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
png_teo_hist = snakemake@output[["png_teo_hist"]] # histogram of teosinte occurence across elevation
# png_teo_hist = "global_ancestry/plots/teosinte_hx_occurence_hist.png"
lm_tex = snakemake@output[["lm_tex"]]
# lm_tex = "global_ancestry/tables/lm_elevation.tex"
png_global_anc_multi = snakemake@output[["png_global_anc_multi"]]
# png_global_anc_multi = "global_ancestry/plots/global_anc_multi.png"

# get colors for plot
# source("colors.R")
source(snakemake@input[["colors"]])

# get NGSAdmix results for K = 2, d_admix2
load(snakemake@input[["k2"]])
# load("global_ancestry/results/NGSAdmix/HILO_MAIZE55/K2_alphas_by_ind.RData")

teosinte_excel = snakemake@input[["teo"]]
# teosinte_excel = "data/zea_occurence/teocintles_historico_jsgetal.xlsx"

# read in teosinte occurence data (parviglumis & mexicana)
teo <- readxl::read_xlsx(teosinte_excel, 
                         sheet = "RegistrosFinal") %>%
  filter(Taxa %in% c("Zea mays parviglumis", "Zea mays mexicana")) %>%
  rename(ELEVATION = Alt,
         zea = Subespecie,
         year = Fecha) %>%
  arrange(ELEVATION) %>%
  mutate(latlonelev = paste(Latitud, Longitud, ELEVATION),
         latlon = paste(Latitud, Longitud),
         latlonelevyear = paste(Latitud, Longitud, ELEVATION, year)) %>%
  filter(!duplicated(latlonelevyear)) %>% # only keep unique occurence observations per location (lat/lon + elevation) and year
  dplyr::select(zea, ELEVATION, Estado, Latitud, Longitud, year)

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
  with(., lm(mexicana ~ ELEVATION))
lm_mex = filter(d_admix2, group == "sympatric_mexicana") %>%
  with(., lm(mexicana ~ ELEVATION))
#glance(lm_maize)
#summary(lm_mex)
#glance(lm_mex)
#summary(lm_maize)

lm_table <- bind_rows(mutate(broom::tidy(lm_mex), subspecies = "maize"),
                      mutate(broom::tidy(lm_mex), subspecies = "mexicana")) %>%
  dplyr::select(subspecies, term, estimate, std.error, statistic, p.value)

# print linear model output to a file
print(xtable(lm_table,
             caption = "\\color{Gray} \\textbf{Elevational ancestry clines} Effect of elevation (m) on genomewide proportion \\texit{mexicana} ancestry in sympatric maize and \\textit{mexicana}",
             label = "anova_lat_clines",
             type = "latex",
             #display = c("f", "s", "s", "g", "g", "g", "g"),
             digits = c(1, 1, 1, -2, -2, 2, -2),
             latex.environments = NULL),
      include.rownames = F,
      file = lm_tex)


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

# plot all individuals (including allopatric) with elevation:
# note: the plotted line still only is fitted to sympatric individuals
shape_group_zea2 <- shape_group_zea 
shape_group_zea2[2] <- 6 # changes to upside down triangle to make more distinguishable
d_admix3 <- d_admix2 %>% # add elevation for Palmar Chico allopatric maize samples
  mutate(., ELEVATION = ifelse(LOCALITY == "Palmar Chico" & 
                                 group == "allopatric_maize", 
                               983, # elevation of Palmar Chico maize landrace samples from Yang 2019 https://doi.org/10.1073/pnas.1820997116
                               ELEVATION)) %>%
  arrange(., ELEVATION) %>%
  mutate(., LOCALITY = ifelse(symp_allo == "allopatric",
                              paste0(LOCALITY, "*"), LOCALITY)) %>% # add astericks to allopatric localities
  mutate(LOCALITY = factor(LOCALITY, ordered = T, levels = unique(.$LOCALITY)))
p_symp_allo_elev <- d_admix3 %>%
  ggplot(., aes(x = ELEVATION, 
                y = mexicana, 
                color = LOCALITY,
                shape = group)) +
  stat_smooth(data = filter(d_admix2, symp_allo == "sympatric"),
              method = "lm", se = T, color = "black") +
  geom_point(data = d_admix3, alpha = 0.75, size = 2) +
  ylab("Proportion mexicana ancestry") +
  xlab("Elevation (m)") +
  ggtitle("Clines in mexicana ancestry across elevation") +
  labs(color = "Location", shape = "Zea") +
  theme_classic() +
  scale_shape_manual(values = shape_group_zea2, labels = zea_group_labels) +
  coord_cartesian(ylim = 0:1, clip = "off") +
  scale_y_continuous(expand=c(0,0)) +
  xlim(c(950, 2650)) + # upper elevation for sympatric pops is 2609m but mexicana grows up to nearly 3000m
  #xlim(c(950, 3000)) +
  #geom_segment(x = 900, xend = 2000, y = -.15, yend = -.15, 
               #color = col_maize_mex_parv[["parviglumis"]],
               #arrow = arrow(ends = "first", type = "open", length = unit(.1, "inches"))) +
  #geom_segment(x = 1500, xend = 3000, y = -.2, yend = -.2, color = col_maize_mex_parv[["mexicana"]]) +
  #theme(plot.margin = unit(c(1,1,5,0), "lines")) +
  scale_color_viridis_d(direction = -1, option = "viridis")
#p_symp_allo_elev
ggsave(filename = png_elev_symp_allo,
       plot = p_symp_allo_elev,
       device = "png", height = 6.5, width = 7.5, units = "in", dpi = 300)

# histogram of teosinte occurence data
p_teo_hist <- teo %>%
  ggplot(.) +
  geom_histogram(aes(x = ELEVATION, fill = zea),
                 alpha = 0.5, 
                 bins = 35,
                 position = "identity") +
  theme_classic() +
  #xlim(c(850, 3050)) +
  scale_fill_manual(values = col_maize_mex_parv) +
  ggtitle("Teosinte occurence data 1842–2016") +
  labs(fill = "Teosinte", x = "Elevation (m)", y = "Observations")
# plot(p_teo_hist)
ggsave(filename = png_teo_hist,
       plot = p_teo_hist,
       device = "png", 
       height = 3, width = 5, 
       units = "in", dpi = 300)

# elevational range limits observed teosinte:
#teo %>%
#  group_by(zea) %>%
#  summarise(min = min(ELEVATION),
#            max = max(ELEVATION))


# multipanel plot of global ancestry data:
# allopatric maize/mexicana structure only
p_structure_allo_maize <- d_admix2 %>%
  filter(group == "allopatric_maize") %>%
  mutate(ELEVATION = ifelse(group == "allopatric_maize", 983, ELEVATION)) %>% # add dummy elevation for allopatric maize (from Palmar Chico)
  arrange(., ELEVATION, popN) %>%
  mutate(sample = 1:nrow(.)) %>%
  tidyr::gather(., "ancestry", "p", c("mexicana", "maize")) %>%
  ggplot(., aes(fill = ancestry, y = p, x = sample)) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(values = col_maize_mex_parv) +
  labs(fill = "Ancestry", y = "Proportion") +
  theme_classic() +
  theme(axis.ticks.x = element_blank()) +
  scale_x_discrete(name = element_blank(), 
                   limits = c(27.5),
                   labels = c("Palmar Chico\n(983m)"))
p_structure_allo_maize
p_structure_allo_mexicana <- d_admix2 %>%
  filter(group == "allopatric_mexicana") %>%
  arrange(., ELEVATION, popN) %>%
  mutate(sample = 1:nrow(.)) %>%
  mutate(sample = ifelse(LOCALITY == "Malinalco", sample + 1, ifelse(LOCALITY == "Amecameca", sample + 2, sample))) %>%
  tidyr::gather(., "ancestry", "p", c("mexicana", "maize")) %>%
  ggplot(., aes(fill = ancestry, y = p, x = sample)) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(values = col_maize_mex_parv) +
  labs(fill = "Ancestry", 
       y = "Proportion") +
  theme_classic() +
  theme(axis.ticks.x = element_blank()) +
  scale_x_discrete(name = element_blank(), 
                   limits = c(5.5, 20, 37.5),
                   labels = c("Puerta Encantada\n(1,658m)", 
                              "Malinalco\n(1,887m)", 
                              "Amecameca\n(2,467m)")) +
  guides(fill = guide_legend(override.aes = list(size = smallPointSize,
                                                 shape = 0.1))) +
  theme(legend.title = element_text(size = smallLegendTitle), 
        legend.text  = element_text(size = smallLegendText),
        legend.key.size = unit(smallLegendSpacing, "lines"),
        legend.key.width = unit(smallLegendSpacing, "lines"))
p_structure_allo_mexicana
# to teo occurrence plot add arrows for our sympatric sampling range
smallPointSize = 0.25
smallLabel = 7
smallLegendText = 6
smallLegendSize = 2
p_combined <- grid.arrange(grobs = list(ggplotGrob(p_symp_elev +
                                       guides(color = F) +
                                       labs(shape = "Sympatric") +
                                       scale_shape_discrete(labels = c("Sympatric maize", "Sympatric mexicana")) +
                                       theme(plot.title = element_blank(),
                                             legend.position = "top",
                                             legend.text = element_text(size = smallLabel + 1),
                                             legend.title = element_blank(),
                                             legend.margin = margin(c(0,0,0,0)))),
                          ggplotGrob(p_structure_allo_mexicana +
                                       labs(subtitle = "Allopatric mexicana") +
                                       guides(fill = F) +
                                       theme(axis.text.x = element_text(size = smallLabel),
                                             plot.subtitle = element_text(size = smallLabel + 1),
                                             axis.title = element_text(size = smallLabel + 1))
                                       ),
                          ggplotGrob(p_structure_allo_maize + 
                                       labs(subtitle = "Allopatric maize") +
                                       guides(fill = F) +
                                       theme(axis.text.x = element_text(size = smallLabel),
                                             plot.subtitle = element_text(size = smallLabel + 1),
                                             axis.title = element_text(size = smallLabel + 1))
                                     ),
                          ggplotGrob(p_teo_hist + #ggtitle("C") + 
                                       theme(plot.title = element_blank(),
                                             axis.title = element_text(size = smallLabel + 1)) +
                                       guides(fill = F) +
                                       geom_segment(aes(x = min(d_admix2$ELEVATION[d_admix2$symp_allo == "sympatric"]), 
                                                        y = 110, 
                                                        xend = max(d_admix2$ELEVATION[d_admix2$symp_allo == "sympatric"]), 
                                                        yend = 110),
                                                    #lineend = "butt",
                                                    #linejoin = "mitre",
                                                    arrow = arrow(length = unit(1, "mm"), 
                                                                  ends = "both",
                                                                  type = "closed")) +
                                       geom_text(aes(x = mean(range(d_admix2$ELEVATION[d_admix2$symp_allo == "sympatric"])), 
                                                     y = 120, 
                                                     label = "sampled range"),
                                                 size = 2)
                                     ),
                          cowplot::get_legend(p_structure_allo_mexicana +
                                                theme(legend.title = element_text(size = smallLabel), 
                                                      legend.text  = element_text(size = smallLegendText),
                                                      legend.key.size = unit(smallLegendSize, "mm"),
                                                      legend.key.width = unit(smallLegendSize, "mm"),
                                                      legend.margin = margin(c(1,1,5,0)),
                                                      plot.margin = margin(c(0,0,0,0)))),
                          cowplot::get_legend(p_teo_hist +
                                                theme(legend.title = element_text(size = smallLabel), 
                                                      legend.text  = element_text(size = smallLegendText),
                                                      legend.key.size = unit(smallLegendSize, "mm"),
                                                      legend.key.width = unit(smallLegendSize, "mm"),
                                                      legend.margin = margin(c(1,1,5,0)),
                                                      plot.margin = margin(c(0,0,0,0)))), # top, right, bottom, left
                          textGrob(label = "A", 
                                   x = unit(0.5, "lines"), 
                                   y = unit(0, "lines")),
                          textGrob(label = "B", 
                                   x = unit(0.5, "lines"), 
                                   y = unit(0, "lines")),
                          textGrob(label = "C", 
                                   x = unit(0.5, "lines"), 
                                   y = unit(0.5, "lines"))),
             layout_matrix = rbind(c(7, 8, NA),
                                   c(1, 2, 5),
                                   c(1, 3, 3),
                                   c(1, 9, NA),
                                   c(1, 4, 6)),
             heights = c(.1, 2, 2, .1, 2),
             widths = c(5, 4.75, 0.75))
#plot_grid(lambda1,lambda12,labels=c('a)','b)'),ncol=1)
p_combined
ggsave(png_global_anc_multi, 
       plot = p_combined, 
       device = "png", 
       width = 7.5, 
       height = 4, 
       units = "in",
       dpi = 300)

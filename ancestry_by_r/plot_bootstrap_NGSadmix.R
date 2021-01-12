#!/usr/bin/env Rscript

# This script plots the the bootstrap results for
# ancestry across recombination quintiles
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(bedr)
library(grid)
library(gridExtra)
library(cowplot)
library(broom)
library(purrr)
library(xtable)

# load variables from Snakefile
alpha = as.numeric(snakemake@params[["alpha"]]) # for confidence intervals
# alpha = 0.05
windows_file = snakemake@input[["windows"]] # feature information for 1cM windows
# windows_file = "ancestry_by_r/results/map_pos_1cM_windows.txt"
inv_file = snakemake@input[["inv"]]
# inv_file = "data/refMaize/inversions/knownInv_v4_coord.txt"

load(snakemake@input[["r5"]]) # NGSAdmix results by recombination quintile
# load("ancestry_by_r/results/bootstrap_1cM/HILO_MAIZE55/r5_K2.Rdata")
load(snakemake@input[["cd5"]]) # NGSAdmix results by coding bp/cM quintile
# load("ancestry_by_r/results/bootstrap_1cM/HILO_MAIZE55/cd5_K2.Rdata")
source(snakemake@input[["colors"]]) # plotting colors
# source("colors.R")
load(snakemake@input[["meta"]]) # meta
# load("samples/HILO_MAIZE55_meta.RData")

# png filenames out:
png_multi = snakemake@output[["png_multi"]]
# png_multi = "ancestry_by_r/plots/K2_by_r_multi_panel.png"
png_r5_symp = snakemake@output[["png_r5_symp"]]
# png_r5_symp = "ancestry_by_r/plots/K2_by_r_bootstrap_sympatric_only.png"
png_r5_symp_allo = snakemake@output[["png_r5_symp_allo"]]
# png_r5_symp_allo = "ancestry_by_r/plots/K2_by_r_bootstrap_sympatric_and_allopatric.png"
png_cd5_symp = snakemake@output[["png_cd5_symp"]]
# png_cd5_symp = "ancestry_by_r/plots/K2_by_cd_bootstrap_sympatric_only.png"
png_cd5_symp_allo = snakemake@output[["png_cd5_symp_allo"]]
# png_cd5_symp_allo = "ancestry_by_r/plots/K2_by_cd_bootstrap_sympatric_and_allopatric.png"
png_cor_r_cd = snakemake@output[["png_cor_r_cd"]]
# png_cor_r_cd = "ancestry_by_r/plots/corr_r_cd_1cM.png"
png_cor_r_frac = snakemake@output[["png_cor_r_frac"]]
# png_cor_r_frac = "ancestry_by_r/plots/corr_r_frac_1cM.png"
png_facet_r5 = snakemake@output[["png_facet_r5"]]
# png_facet_r5 = "ancestry_by_r/plots/K2_by_r_bootstrap_lm_elevation_facet_r.png"
png_color_elev_r5 = snakemake@output[["png_color_elev_r5"]]
# png_color_elev_r5 = "ancestry_by_r/plots/K2_by_r_bootstrap_lm_elevation_color_elev.png"
file_elev_r_interaction = snakemake@output[["file_elev_r_interaction"]]
# file_elev_r_interaction = "ancestry_by_r/tables/elev_r_interaction.tex"

# load inversion coordinates
inv = read.table(inv_file, stringsAsFactors = F, header = F) %>%
  data.table::setnames(c("ID", "chr", "start", "end", "length"))

# load windows
windows <- read.table(windows_file, header = T, stringsAsFactors = F, sep = "\t") %>%
  dplyr::arrange(quintile_r5) %>%
  dplyr::mutate(bin_r5 = factor(bin_r5, levels = unique(bin_r5), ordered = T)) %>%
  dplyr::arrange(quintile_cd5) %>%
  dplyr::mutate(bin_cd5 = factor(bin_cd5, levels = unique(bin_cd5), ordered = T)) %>%
  dplyr::arrange(quintile_frac5) %>%
  dplyr::mutate(bin_frac5 = factor(bin_frac5, levels = unique(bin_frac5), ordered = T)) %>%
  dplyr::arrange(chr, start)

windows_by_inv = bedr(engine = "bedtools",
     input = list(a = windows %>%
                    dplyr::select("chr", "start", "end", "window") %>%
                    dplyr::mutate(., chr = as.character(chr)), 
                  b = dplyr::mutate(inv, chr = as.character(chr)) %>%
                    dplyr::select(chr, start, end, ID) %>%
                    dplyr::arrange(chr, start)),
     method = "map",
     params = "-o distinct -c 4",
     check.chr = F) %>%
  data.table::setnames(c("chr", "start", "end", "window", "inv")) %>%
  dplyr::select(window, inv) %>%
  dplyr::mutate(inv = ifelse(inv == ".", "None", inv))


# metadata by population for local ancestry
meta_symp = dplyr::filter(meta, symp_allo == "sympatric") %>%
  dplyr::filter(., est_coverage >= 0.5) %>%
  dplyr::group_by(popN, RI_ACCESSION, zea, symp_allo, group, ELEVATION, LOCALITY, GEOCTY) %>%
  dplyr::summarise(n_local_ancestry = n()) %>%
  dplyr::mutate(pop = paste0("pop", popN))

# make plots
# mexicana ancestry for K=2 in sympatric maize and mexicana:
p_r5_symp <- r5$anc_boot_perc %>%
  dplyr::filter(ancestry == "mexicana_ancestry") %>%
  dplyr::filter(symp_allo == "sympatric") %>%
  ggplot(aes(x = bin, y = p, group = zea)) +
  # first plot original point estimates for ind. ancestry
  geom_point(data = filter(r5$anc_ind, ancestry == "mexicana"),
             aes(x = bin,
                 y = p,
                 shape = zea,
                 color = zea),
             position = position_jitter(0.2)) +
  scale_color_manual(values = col_maize_mex_parv) +
  # then add mean for that group
  geom_point(pch = 18, size = 2) +
  # and errorbars for 95% CI around that mean
  # based on bootstrap with NGSAdmix (percentile method)
  geom_errorbar(aes(ymin = low_boot,
                    ymax = high_boot),
                width = .5) +
  xlab("Recombination rate quintile (cM/Mb)") +
  ylab("Proportion mexicana ancestry") +
  theme_classic() +
  guides(color = guide_legend("Subspecies"),
         shape = guide_legend("Subspecies")) +
  ggtitle("Mexicana ancestry by recombination rate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p_r5_symp
ggsave(file = png_r5_symp,
       plot = p_r5_symp,
       device = "png",
       width = 5, height = 4, 
       units = "in", dpi = 300)

# by coding bp per cM
p_cd5_symp <- cd5$anc_boot_perc %>%
  dplyr::filter(ancestry == "mexicana_ancestry") %>%
  dplyr::filter(symp_allo == "sympatric") %>%
  ggplot(aes(x = bin, y = p, group = zea)) +
  # first plot original point estimates for ind. ancestry
  geom_point(data = filter(cd5$anc_ind, zea != "parviglumis" & ancestry == "mexicana"),
             aes(x = bin,
                 y = p,
                 shape = zea,
                 color = zea),
             position = position_jitter(0.2)) +
  scale_color_manual(values = col_maize_mex_parv) +
  # then add mean for that group
  geom_point(pch = 18, size = 2) +
  # and errorbars for 95% CI around that mean
  # based on bootstrap with NGSAdmix (based on percentile method)
  geom_errorbar(aes(ymin = low_boot,
                    ymax = high_boot),
                width = .5) +
  xlab("Gene density (coding bp/cM)") +
  ylab("Proportion mexicana ancestry") +
  theme_classic() +
  guides(color = guide_legend("Subspecies"),
         shape = guide_legend("Subspecies")) +
  ggtitle("Mexicana ancestry by gene density") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p_cd5_symp
ggsave(file = png_cd5_symp,
       plot = p_cd5_symp,
       device = "png",
       width = 5, height = 4, 
       units = "in", dpi = 300)

# mexicana ancestry for K=2 in sympatric and allopatric maize/mex
p_r5_symp_allo <- r5$anc_boot_perc %>%
  dplyr::filter(ancestry == "mexicana_ancestry") %>%
  ggplot(aes(x = bin, y = p, group = zea)) +
  # first plot original point estimates for ind. ancestry
  geom_point(data = filter(r5$anc_ind, zea != "parviglumis" & ancestry == "mexicana"),
             aes(x = bin,
                 y = p,
                 shape = zea,
                 color = zea),
             position = position_jitter(0.2)) +
  scale_color_manual(values = col_maize_mex_parv) +
  # then add mean for that group
  geom_point(pch = 18, size = 2) +
  # and errorbars for 95% CI around that mean
  # based on bootstrap with NGSAdmix (percentile method)
  geom_errorbar(aes(ymin = low_boot,
                    ymax = high_boot),
                width = .5) +
  facet_wrap(~symp_allo) +
  theme_classic() +
  xlab("Recombination rate quintile (cM/Mb)") +
  ylab("Proportion mexicana ancestry") +
  guides(color = guide_legend("Subspecies", override.aes = list(size = 3)),
         shape = guide_legend("Subspecies")) +
  #ggtitle("Mexicana ancestry by recombination rate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p_r5_symp_allo
ggsave(file = png_r5_symp_allo,
       plot = p_r5_symp_allo,
       device = "png",
       width = 7, height = 7, 
       units = "in", dpi = 300)

p_cd5_symp_allo <- cd5$anc_boot_perc %>%
  dplyr::filter(ancestry == "mexicana_ancestry") %>%
  ggplot(aes(x = bin, y = p, group = zea)) +
  # first plot original point estimates for ind. ancestry
  geom_point(data = filter(cd5$anc_ind, zea != "parviglumis" & ancestry == "mexicana"),
             aes(x = bin,
                 y = p,
                 shape = zea,
                 color = zea),
             position = position_jitter(0.2)) +
  scale_color_manual(values = col_maize_mex_parv) +
  # then add mean for that group
  geom_point(pch = 18, size = 2) +
  # and errorbars for 95% CI around that mean
  # based on bootstrap with NGSAdmix (percentile method)
  geom_errorbar(aes(ymin = low_boot,
                    ymax = high_boot),
                width = .5) +
  xlab("Gene density (coding bp/cM)") +
  ylab("Proportion mexicana ancestry") +
  theme_classic() +
  facet_wrap(~symp_allo) +
  guides(color = guide_legend("Subspecies", override.aes = list(size = 3)),
         shape = guide_legend("Subspecies")) +
  #ggtitle("Mexicana ancestry by gene density") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p_cd5_symp_allo
ggsave(file = png_cd5_symp_allo,
       plot = p_cd5_symp_allo,
       device = "png",
       width = 7, height = 7, 
       units = "in", dpi = 300)


# plot feature correlations for genomic windows
p_cor_r_frac <- ggplot(windows, aes(x = log10(cM_Mb), y = frac_bp_coding, 
                    col = bin_frac5, shape = bin_r5)) +
  geom_point() +
  ggtitle("Log(cM/Mb) vs. fraction coding (vs. noncoding) bp in 1cM windows")
#p_cor_r_frac
ggsave(file = png_cor_r_frac,
       plot = p_cor_r_frac,
       device = "png",
       width = 7, height = 4, 
       units = "in", dpi = 300)
print("Genomic feature correlation at 1cM window scale:")
print("cM_Mb ~ frac_bp_coding")
print(cor(windows$cM_Mb, windows$frac_bp_coding))
print("quintile_r5 ~ quintile_frac5")
print(cor(windows$quintile_r5, windows$quintile_frac5))
print("log10(cM_Mb) ~ frac_bp_coding")
print(cor(log10(windows$cM_Mb), windows$frac_bp_coding))

p_cor_r_cd <- ggplot(windows, aes(x = log10(cM_Mb), y = coding_bp, 
                    col = bin_cd5, shape = bin_r5)) +
  geom_point() +
  #ggtitle("Log(cM/Mb) vs. coding bp in 1cM windows") +
  labs(color = "Gene density quintile\n(coding bp/cM)",
       shape = "Recombination rate quintile\n(cM/Mb)",
       y = "Coding bp/cM",
       x = "log10(cM/Mb)")
#p_cor_r_cd
ggsave(file = png_cor_r_cd,
       plot = p_cor_r_cd,
       device = "png",
       width = 6, height = 4, 
       units = "in", dpi = 300)
print("log10(cM_Mb_ ~ coding_bp")
print(cor(log10(windows$cM_Mb), windows$coding_bp))
print("quintile_r5 ~ quintile_cd5")
print(cor(windows$quintile_r5, windows$quintile_cd5))

# plot ancestry across recombination bins and elevation:
p_facet_r5 <- r5$anc_ind %>%
  dplyr::filter(., symp_allo == "sympatric") %>%
  dplyr::filter(., ancestry == "mexicana") %>%
  ggplot(., aes(x = ELEVATION, y = p,
                color = zea,
                shape = zea)) +
  geom_point() +
  geom_smooth(method = "lm", aes(group = zea, color = zea, fill = zea)) +
  #ggtitle("Mexicana ancestry by recombination rate and elevation") +
  scale_color_manual(values = col_maize_mex_parv) +
  scale_fill_manual(values = col_maize_mex_parv) +
  theme_light() +
  facet_wrap(~ bin) +
  xlab("Elevation (m)") +
  ylab("Proportion mexicana ancestry") +
  guides(color = guide_legend("Subspecies"),
         shape = guide_legend("Subspecies"),
         fill = guide_legend("Subspecies"))
#p_facet_r5
ggsave(file = png_facet_r5,
       plot = p_facet_r5,
       device = "png",
       width = 7, height = 4, 
       units = "in", dpi = 300)


# different way of visualizing the same pattern:
p_color_elev_r5 <- r5$anc_ind %>%
  dplyr::filter(., symp_allo == "sympatric") %>%
  dplyr::filter(., ancestry == "mexicana") %>%
  ggplot(., aes(x = bin, y = p,
                color = ELEVATION)) +
  geom_smooth(method = "lm", se = FALSE, 
              aes(group = popN, color = ELEVATION)) +
  #ggtitle("Mexicana ancestry by recombination rate and elevation") +
  scale_color_viridis_c(direction = -1, option = "viridis") +
  theme_light() +
  facet_wrap(~zea) +
  ylab("Proportion mexicana ancestry") +
  xlab("Recombination rate quintile (cM/Mb)") +
  guides(color = guide_legend("Elevation (m)")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p_color_elev_r5
ggsave(file = png_color_elev_r5,
       plot = p_color_elev_r5,
       device = "png",
       width = 7, height = 6, 
       units = "in", dpi = 300)

# multipanel plot:
# facet label names
r5.labs <- c("Low r: [0.0185, 0.0524] cM/Mb", 
             "High r: (1.09, 287] cM/Mb")
names(r5.labs) <- c("[0.0185,0.0524]", "(1.09,287]")

# plots across elevation of lowest and highest quintiles
p_elev_facet_2 <- r5$anc_ind %>%
  dplyr::filter(., symp_allo == "sympatric") %>%
  dplyr::filter(., ancestry == "mexicana") %>%
  dplyr::filter(., quintile %in% c(1,5)) %>% # only show lowest and highest recomb. quintiles
  ggplot(., aes(x = ELEVATION, y = p,
                color = zea,
                shape = zea)) +
  geom_point() +
  geom_smooth(method = "lm", aes(group = zea, color = zea, fill = zea)) +
  scale_color_manual(values = col_maize_mex_parv) +
  scale_fill_manual(values = col_maize_mex_parv) +
  theme_light() +
  facet_wrap(~bin, nrow = 2, labeller = labeller(bin = r5.labs)) +
  xlab("Elevation (m)") +
  ylab("Proportion mexicana ancestry") +
  guides(color = guide_legend("Subspecies"),
         shape = guide_legend("Subspecies"),
         fill = guide_legend("Subspecies"))
#p_elev_facet_2

p_multi <- grid.arrange(grobs = list(ggplotGrob(p_r5_symp_allo +
                                                  theme(legend.position = "none")),
                                        ggplotGrob(p_elev_facet_2 +
                                                     theme(legend.position = "none",
                                                           plot.margin = margin(l = 8, r = 5.5, 
                                                                                t = 5.5, b = 5.5, 
                                                                                unit = "pt"))),
                                     textGrob(label = "A", 
                                              x = unit(0.5, "lines"), 
                                              y = unit(0, "lines")),
                                     textGrob(label = "B", 
                                              x = unit(0.5, "lines"), 
                                              y = unit(0, "lines")),
                                     cowplot::get_legend(p_r5_symp_allo +
                                                           theme(legend.position = "bottom"))
                                        ),
                           layout_matrix = rbind(c(3,4), c(1,2), c(5,5)),
                           heights = c(0.1, 1, 0.1),
                           widths = c(5, 3))
p_multi
ggsave(png_multi, 
       plot = p_multi, 
       device = "png", 
       width = 7.5, 
       height = 6, 
       units = "in",
       dpi = 300)

# test for significant difference in slope ancestry ~ elev
# between highest and lowest r bins
tbl_elev_r_interation = r5$anc_ind %>%
  dplyr::filter(., symp_allo == "sympatric") %>%
  dplyr::filter(., ancestry == "mexicana") %>%
  dplyr::filter(., quintile %in% c(1,5)) %>%
  mutate(.,
    elevation_km = ELEVATION/1000,
    r = ifelse(quintile == 5, "high_r", ifelse(quintile == 1, "low_r", NA))) %>%
  nest(., -zea) %>%
  mutate(.,
         model = map(data, ~lm(p ~ elevation_km + r + r*elevation_km, 
                               data = .)),
         tidy = map(model, tidy)) %>%
  unnest(tidy) %>%
  dplyr::select(-model, -data) %>%
  mutate(model = "mexicana ancestry ~ elevation + r + r*elevation") %>%
  rename(term_messy = term) %>%
  # make a pretty table
  left_join(., data.frame(term_messy = c("(Intercept)", "elevation_km",
                                         "rlow_r", "elevation_km:rlow_r"),
                          term = c("intercept", "elevation (km)",
                                   "r (low recombination)", "elevation*r"),
                          stringsAsFactors = F),
            by = "term_messy") %>%
  dplyr::select(zea, term, estimate, std.error, statistic, p.value) 

# print table to file
print(xtable(tbl_elev_r_interation, 
             digits = c(1, 0, 0, 3, 3, 3, -2),
             label = "tbl_elev_r_interaction",
             type = "latex", 
             latex.environments = NULL),
      include.rownames = F,
      file = file_elev_r_interaction)
      
  
# use all quintiles (ordinal) as numeric
# and interaction term is also significant
print("maize then mexicana model results -- using ordinal r quintiles as numeric")
for (z in c("maize", "mexicana")){
  r5$anc_ind %>%
    dplyr::filter(., symp_allo == "sympatric") %>%
    dplyr::filter(., ancestry == "mexicana") %>%
    mutate(., 
           elevation_km = ELEVATION/1000) %>%
    dplyr::filter(., zea == z) %>%
    with(data = ., lm(p ~ elevation_km + quintile + quintile*elevation_km)) %>%
    summary(.) %>%
    print(.)
}

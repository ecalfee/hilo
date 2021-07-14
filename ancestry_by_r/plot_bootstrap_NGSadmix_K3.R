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
library(stringr)

# load variables from Snakefile
windows_file = snakemake@input[["windows"]] # feature information for 1cM windows
# windows_file = "ancestry_by_r/results/map_pos_1cM_windows.txt"
inv_file = snakemake@input[["inv"]]
# inv_file = "data/refMaize/inversions/knownInv_v4_coord.txt"
alpha = 0.05 # use 95% confidence level for bootstraps
load(snakemake@input[["r5"]]) # NGSAdmix results by recombination quintile
# load("ancestry_by_r/results/bootstrap_1cM/HILO_MAIZE55_PARV50/r5_K3.Rdata")
load(snakemake@input[["cd5"]]) # NGSAdmix results by coding bp/cM quintile
# load("ancestry_by_r/results/bootstrap_1cM/HILO_MAIZE55_PARV50/cd5_K3.Rdata")
source(snakemake@input[["colors"]]) # plotting colors
# source("colors.R")
load(snakemake@input[["meta"]]) # meta
# load("samples/HILO_MAIZE55_PARV50_meta.RData")

# png filenames out:
png_multi = snakemake@output[["png_multi"]]
# png_multi = "ancestry_by_r/plots/HILO_MAIZE55_PARV50_K3_by_r_multi_panel.png"
png_multi_lzw = snakemake@output[["png_multi_lzw"]]
# png_multi_lzw = "../hilo_manuscript/figures_main/HILO_MAIZE55_PARV50_K3_by_r_multi_panel.tif"

png_r5_symp = snakemake@output[["png_r5_symp"]]
# png_r5_symp = "ancestry_by_r/plots/HILO_MAIZE55_PARV50_K3_by_r_bootstrap_sympatric_only.png"
png_r5_symp_allo = snakemake@output[["png_r5_symp_allo"]]
# png_r5_symp_allo = "ancestry_by_r/plots/HILO_MAIZE55_PARV50_K3_by_r_bootstrap_sympatric_and_allopatric.png"
png_r5_symp_allo_lzw = snakemake@output[["png_r5_symp_allo_lzw"]]
# png_r5_symp_allo_lzw = "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_by_r_bootstrap_sympatric_and_allopatric.tif"

png_cd5_symp = snakemake@output[["png_cd5_symp"]]
# png_cd5_symp = "ancestry_by_r/plots/HILO_MAIZE55_PARV50_K3_by_cd_bootstrap_sympatric_only.png"
png_cd5_symp_allo = snakemake@output[["png_cd5_symp_allo"]]
# png_cd5_symp_allo = "ancestry_by_r/plots/HILO_MAIZE55_PARV50_K3_by_cd_bootstrap_sympatric_and_allopatric.png"
png_cd5_symp_allo_lzw = snakemake@output[["png_cd5_symp_allo_lzw"]]
# png_cd5_symp_allo_lzw = "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_by_cd_bootstrap_sympatric_and_allopatric.tif"

png_cor_r_cd = snakemake@output[["png_cor_r_cd"]]
# png_cor_r_cd = "ancestry_by_r/plots/HILO_MAIZE55_PARV50_K3_corr_r_cd_1cM.png"
png_cor_r_frac = snakemake@output[["png_cor_r_frac"]]
# png_cor_r_frac = "ancestry_by_r/plots/HILO_MAIZE55_PARV50_K3_corr_r_frac_1cM.png"

png_color_elev_r5 = snakemake@output[["png_color_elev_r5"]]
# png_color_elev_r5 = "ancestry_by_r/plots/HILO_MAIZE55_PARV50_K3_by_r_bootstrap_lm_elevation_color_elev.png"
png_color_elev_r5_lzw = snakemake@output[["png_color_elev_r5_lzw"]]
# png_color_elev_r5_lzw = "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_by_r_bootstrap_lm_elevation_color_elev.tif"

file_elev_r_interaction = snakemake@output[["file_elev_r_interaction"]] # high vs. low
# file_elev_r_interaction = "ancestry_by_r/tables/HILO_MAIZE55_PARV50_K3_elev_r_interaction.tex"

file_elev_r_interaction_5 = snakemake@output[["file_elev_r_interaction_5"]] # all 5 r bins (for supplement)
# file_elev_r_interaction_5 = "ancestry_by_r/tables/HILO_MAIZE55_PARV50_K3_elev_r_interaction_5.tex"
file_elev_r_interaction_5_tex = snakemake@output[["file_elev_r_interaction_5_tex"]] # all 5 r bins (for supplement)
# file_elev_r_interaction_5_tex = "../hilo_manuscript/tables/HILO_MAIZE55_PARV50_K3_elev_r_interaction_5.tex"

file_spearmans_rho_ngsadmix = snakemake@output[["file_spearmans_rho_ngsadmix"]]
# file_spearmans_rho_ngsadmix = "ancestry_by_r/tables/HILO_MAIZE55_PARV50_K3_spearmans_rho_ngsadmix.tex"
file_spearmans_rho_ngsadmix_tex = snakemake@output[["file_spearmans_rho_ngsadmix_tex"]]
# file_spearmans_rho_ngsadmix_tex = "../hilo_manuscript/tables/HILO_MAIZE55_PARV50_K3_spearmans_rho_ngsadmix.tex"

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
# K=3 introgression ~ r
p_r5_symp = r5$anc_boot_perc %>%
  dplyr::filter(symp_allo == "sympatric") %>%
  dplyr::mutate(ancestry = stringr::str_extract(ancestry, "[a-z]+")) %>% # makes 'maize_ancestry' into just 'maize'
  dplyr::filter(ancestry != zea) %>% # don't plot e.g. maize ancestry within maize
  dplyr::mutate(label_introgression = factor(paste(ancestry, "in\n", symp_allo, zea),
                                             ordered = T,
                levels = paste(c("mexicana", "maize", "parviglumis", "parviglumis"),
                               "in\n", "sympatric",
                               c("maize", "mexicana", "maize", "mexicana")))) %>%
  ggplot(aes(x = bin, y = p, group = zea)) +
  # first plot original point estimates for ind. ancestry
  geom_point(data = filter(r5$anc_ind,
                           symp_allo == "sympatric" & ancestry != zea) %>%
               dplyr::mutate(label_introgression = factor(paste(ancestry, "in\n", symp_allo, zea),
                                                          ordered = T)),
             aes(x = bin,
                 y = p,
                 shape = zea,
                 color = ancestry),
             alpha = 0.6,
             size = 1,
             position = position_jitter(0.2)) +
  scale_color_manual(values = col_maize_mex_parv) +
  # then add mean for that group
  geom_point(pch = 18, size = 2) +
  # and errorbars for 95% CI around that mean
  # based on bootstrap with NGSAdmix (percentile method)
  geom_errorbar(aes(ymin = low_boot,
                    ymax = high_boot),
                width = .5) +
  facet_wrap(~ label_introgression) +
  theme_classic() +
  xlab("Recombination rate quintile (cM/Mb)") +
  ylab("Introgressed ancestry proportion") +
  guides(color = guide_legend("Ancestry", override.aes = list(size = 3)),
         shape = guide_legend("Ancestry")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None")
# p_r5_symp
ggsave(file = png_r5_symp,
       plot = p_r5_symp,
       device = "png",
       width = 5, height = 4,
       units = "in", dpi = 300)

# plot sympatric and reference together: NGSAdmix ancestry ~ r
p_r5_symp_allo = r5$anc_boot_perc %>%
  dplyr::mutate(ref_symp = ifelse(zea == "parviglumis" | symp_allo == "allopatric",
                "reference populations", "sympatric populations")) %>%
  #dplyr::filter(symp_allo == "sympatric") %>%
  dplyr::mutate(ancestry = stringr::str_extract(ancestry, "[a-z]+")) %>% # makes 'maize_ancestry' into just 'maize'
  dplyr::filter(zea != ancestry) %>% # introgressed ancestry only
  ggplot(aes(x = bin, y = p, group = zea, color = ancestry, shape = zea)) +
  # plot mean per group
  geom_point() +
  # and errorbars for 95% CI around that mean
  # based on bootstrap with NGSAdmix (percentile method)
  geom_errorbar(aes(ymin = low_boot,
                    ymax = high_boot),
                width = 0.5) +
  scale_color_manual(values = col_maize_mex_parv) +
  #scale_shape_manual(values = c(19,17,0)) +
  facet_grid(ref_symp ~ zea, scales="free_y") +
  theme_light() +
  xlab("Recombination rate quintile (cM/Mb)") +
  ylab("Introgressed ancestry proportion") +
  guides(color = guide_legend("Introgressed ancestry",
                              override.aes = list(size = 2, linetype = 0)),
         shape = guide_legend("Population subspecies"),
         fill = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p_r5_symp_allo
ggsave(file = png_r5_symp_allo,
       plot = p_r5_symp_allo,
       device = "png",
       width = 7, height = 7,
       units = "in", dpi = 300)
ggsave(file = png_r5_symp_allo_lzw,
       plot = p_r5_symp_allo,
       device = "tiff",
       width = 7, height = 7,
       units = "in", dpi = 300,
       compression = "lzw", type = "cairo")

# by coding bp per cM
p_cd5_symp = cd5$anc_boot_perc %>%
  dplyr::filter(symp_allo == "sympatric") %>%
  dplyr::mutate(ancestry = stringr::str_extract(ancestry, "[a-z]+")) %>% # makes 'maize_ancestry' into just 'maize'
  dplyr::filter(ancestry != zea) %>% # don't plot e.g. maize ancestry within maize
  dplyr::mutate(label_introgression = factor(paste(ancestry, "in\n", symp_allo, zea),
                                             ordered = T,
                                             levels = paste(c("mexicana", "maize", "parviglumis", "parviglumis"),
                                                            "in\n", "sympatric",
                                                            c("maize", "mexicana", "maize", "mexicana")))) %>%
  ggplot(aes(x = bin, y = p, group = zea)) +
  # first plot original point estimates for ind. ancestry
  geom_point(data = filter(cd5$anc_ind,
                           symp_allo == "sympatric" & ancestry != zea) %>%
               dplyr::mutate(label_introgression = factor(paste(ancestry, "in\n", symp_allo, zea),
                                                          ordered = T)),
             aes(x = bin,
                 y = p,
                 shape = zea,
                 color = ancestry),
             alpha = 0.6,
             size = 1,
             position = position_jitter(0.2)) +
  scale_color_manual(values = col_maize_mex_parv) +
  # then add mean for that group
  geom_point(pch = 18, size = 2) +
  # and errorbars for 95% CI around that mean
  # based on bootstrap with NGSAdmix (percentile method)
  geom_errorbar(aes(ymin = low_boot,
                    ymax = high_boot),
                width = .5) +
  facet_wrap(~ label_introgression) +
  theme_classic() +
  xlab("Gene density (coding bp/cM)") +
  ylab("Introgressed ancestry proportion") +
  guides(color = guide_legend("Ancestry", override.aes = list(size = 3)),
         shape = guide_legend("Ancestry")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None")
# p_cd5_symp
ggsave(file = png_cd5_symp,
       plot = p_cd5_symp,
       device = "png",
       width = 5, height = 4,
       units = "in", dpi = 300)


# sympatric and reference populations together: NGSAdmix ancestry ~ coding density
p_cd5_symp_allo = cd5$anc_boot_perc %>%
  dplyr::mutate(ref_symp = ifelse(zea == "parviglumis" | symp_allo == "allopatric",
                                  "reference populations", "sympatric populations")) %>%
  dplyr::mutate(ancestry = stringr::str_extract(ancestry, "[a-z]+")) %>% # makes 'maize_ancestry' into just 'maize'
  dplyr::filter(zea != ancestry) %>% # introgressed ancestry only
  ggplot(aes(x = bin, y = p, group = zea, color = ancestry, shape = zea)) +
  # plot mean per group
  geom_point() +
  # and errorbars for 95% CI around that mean
  # based on bootstrap with NGSAdmix (percentile method)
  geom_errorbar(aes(ymin = low_boot,
                    ymax = high_boot),
                width = 0.5) +
  scale_color_manual(values = col_maize_mex_parv) +
  #scale_shape_manual(values = c(19,17,0)) +
  facet_grid(ref_symp ~ zea, scales="free_y") +
  theme_light() +
  xlab("Recombination rate quintile (cM/Mb)") +
  ylab("Gene density (coding bp/cM)") +
  guides(color = guide_legend("Introgressed ancestry",
                              override.aes = list(size = 2, linetype = 0)),
         shape = guide_legend("Population subspecies"),
         fill = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p_cd5_symp_allo
ggsave(file = png_cd5_symp_allo,
       plot = p_cd5_symp_allo,
       device = "png",
       width = 7, height = 7,
       units = "in", dpi = 300)
ggsave(file = png_cd5_symp_allo_lzw,
       plot = p_cd5_symp_allo,
       device = "tiff",
       width = 7, height = 7,
       units = "in", dpi = 300,
       compression = "lzw", type = "cairo")

# plot feature correlations for genomic windows
p_cor_r_frac <- ggplot(windows, aes(x = log10(cM_Mb), y = frac_bp_coding,
                    col = bin_frac5, shape = bin_r5)) +
  geom_point() +
  ggtitle("Log(cM/Mb) vs. fraction coding (vs. noncoding) bp in 1cM windows")
# p_cor_r_frac
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
# p_cor_r_cd
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
# for sympatric maize
p_elev_facet_r5_maize <- r5$anc_ind %>%
  dplyr::filter(., symp_allo == "sympatric") %>%
  dplyr::filter(., ancestry == "mexicana") %>%
  dplyr::filter(., zea == "maize") %>%
  ggplot(., aes(x = ELEVATION, y = p,
                color = ancestry,
                shape = zea)) +
  geom_point() +
  geom_smooth(method = "lm", aes(group = zea, color = ancestry, fill = ancestry)) +
  #ggtitle("Mexicana ancestry by recombination rate and elevation") +
  scale_color_manual(values = col_maize_mex_parv) +
  scale_fill_manual(values = col_maize_mex_parv) +
  theme_light() +
  facet_wrap(~ bin, ncol = 5) +
  xlab("Elevation (m)") +
  ylab("Mexicana ancestry proportion") +
  guides(color = guide_legend("Subspecies"),
         shape = guide_legend("Subspecies"),
         fill = guide_legend("Subspecies")) +
  theme(legend.position = "None")
# p_elev_facet_r5_maize

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
# p_color_elev_r5
ggsave(file = png_color_elev_r5,
       plot = p_color_elev_r5,
       device = "png",
       width = 7, height = 6,
       units = "in", dpi = 300)
ggsave(file = png_color_elev_r5_lzw,
       plot = p_color_elev_r5,
       device = "tiff",
       width = 7, height = 6,
       units = "in", dpi = 300,
       compression = "lzw", type = "cairo")

# multipanel plot:
# facet label names
r5.labs0 <- paste0(c("lowest r", "low r", "middle r", "high r", "highest r"))
r5.labs <- paste(levels(r5$anc_ind$bin), "cM/Mb")
names(r5.labs) <- levels(r5$anc_ind$bin)
p_r5_symp_maize_mex = p_r5_symp = r5$anc_boot_perc %>%
  dplyr::filter(symp_allo == "sympatric") %>%
  dplyr::mutate(ancestry = stringr::str_extract(ancestry, "[a-z]+")) %>% # makes 'maize_ancestry' into just 'maize'
  dplyr::filter(ancestry != zea & ancestry != "parviglumis") %>% # don't plot e.g. maize ancestry within maize
  dplyr::mutate(label_introgression = factor(paste(ancestry, "ancestry in\n", symp_allo, zea),
                                             ordered = T,
                                             levels = paste(c("mexicana", "maize", "parviglumis", "parviglumis"),
                                                            "ancestry in\n", "sympatric",
                                                            c("maize", "mexicana", "maize", "mexicana")))) %>%
  ggplot(aes(x = bin, y = p, group = zea)) +
  # first plot original point estimates for ind. ancestry
  geom_point(data = filter(r5$anc_ind,
                           symp_allo == "sympatric" & ancestry != zea & ancestry != "parviglumis") %>%
               dplyr::mutate(label_introgression = factor(paste(ancestry, "ancestry in\n", symp_allo, zea),
                                                          ordered = T)),
             aes(x = bin,
                 y = p,
                 shape = zea,
                 color = ancestry),
             alpha = 0.6,
             size = 1,
             position = position_jitter(0.2)) +
  scale_color_manual(values = col_maize_mex_parv) +
  # then add mean for that group
  geom_point(pch = 18, size = 2) +
  # and errorbars for 95% CI around that mean
  # based on bootstrap with NGSAdmix (percentile method)
  geom_errorbar(aes(ymin = low_boot,
                    ymax = high_boot),
                width = .5) +
  facet_wrap(~ label_introgression) +
  theme_classic() +
  xlab("Recombination rate quintile (cM/Mb)") +
  ylab("Introgressed ancestry proportion") +
  guides(color = guide_legend("Ancestry", override.aes = list(size = 3)),
         shape = guide_legend("Ancestry")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None")
p_r5_symp_maize_mex




p_multi <- grid.arrange(grobs = list(textGrob(label = "A",
                                              just = "top",
                                              x = unit(0.5, "lines")),
                                     ggplotGrob(p_r5_symp_maize_mex),
                                     textGrob(label = "B",
                                              just = "top",
                                              x = unit(0.5, "lines"),
                                              y = unit(1.5, "lines")),
                                     ggplotGrob(p_elev_facet_r5_maize +
                                                  theme(legend.position = "none") +
                                                  facet_wrap(~ bin, ncol = 5, labeller = labeller(bin = r5.labs)))

                                     ),
                           #layout_matrix = rbind(c(3, NA, 4, NA),
                          #                       c(NA, 1, NA, 2)),
                           nrow = 4,
                           heights = c(.3, 5, .3, 2.5))


#p_multi
ggsave(png_multi,
       plot = p_multi,
       device = "png",
       width = 7.5,
       height = 6.5,
       units = "in",
       dpi = 300)
ggsave(png_multi_lzw,
       plot = p_multi,
       device = "tiff",
       width = 7.5,
       height = 6.5,
       units = "in",
       dpi = 300,
       compression = "lzw", type = "cairo")

# test for significant difference in slope ancestry ~ elev
# using all quintiles as a numeric scale 0-4
tbl_elev_r_interaction_5 = r5$anc_ind %>%
  dplyr::filter(., symp_allo == "sympatric") %>%
  dplyr::filter(., ancestry == "mexicana") %>%
  mutate(.,
         quintile = quintile - 1, # make quintiles 0-4 instead of 1-5
         elevation_km = ELEVATION/1000) %>%
  nest(., -zea) %>%
  mutate(.,
         model = map(data,
                     ~lm(p ~ elevation_km + quintile + quintile*elevation_km,
                               data = .)),
         tidy = map(model, tidy)) %>%
  unnest(tidy) %>%
  dplyr::select(-model, -data) %>%
  mutate(model = "mexicana ancestry ~ elevation + quintile + quintile*elevation") %>%
  rename(term_messy = term) %>%
  # make a pretty table
  left_join(., data.frame(term_messy = c("(Intercept)", "elevation_km",
                                         "quintile", "elevation_km:quintile"),
                          term = c("intercept", "elevation (km)",
                                   "r quintile", "elevation*r quintile"),
                          stringsAsFactors = F),
            by = "term_messy") %>%
  dplyr::mutate(group = paste("sympatric", zea)) %>%
  dplyr::select(group, term, estimate, std.error, statistic, p.value)

# print table to file
print(xtable(tbl_elev_r_interaction_5,
             digits = c(1, 0, 0, 3, 3, 3, -2),
             label = "tbl_elev_r_interaction_5",
             type = "latex",
             latex.environments = NULL),
      include.rownames = F,
      file = file_elev_r_interaction_5)
print(xtable(tbl_elev_r_interaction_5,
             digits = c(1, 0, 0, 3, 3, 3, -2),
             label = "tbl_elev_r_interaction_5",
             type = "latex",
             latex.environments = NULL),
      include.rownames = F,
      file = file_elev_r_interaction_5_tex)

# make tidy table of rho spearman's rank correlations across quintiles
rho = bind_rows(mutate(r5$spearman,
                       method = "NGSAdmix",
                       feature = "recombination rate (cM/Mb)"),
                mutate(cd5$spearman,
                       method = "NGSAdmix",
                       feature = "gene density (coding bp/cM)")) %>%
  mutate(group = stringr::str_replace(group, "_", " "),
         resolution = "genomic quintiles") %>%
  dplyr::rename(label = group) %>%
  left_join(.,
            data.frame(group = c("sympatric maize", "sympatric mexicana",
                                 "reference parviglumis",
                                 "reference maize", "reference mexicana"),
                       label = c("sympatric maize", "sympatric mexicana",
                                 "parviglumis",
                                 "allopatric maize", "allopatric mexicana"),
                       stringsAsFactors = F),
            by = "label") %>%
  dplyr::select(group, feature, ancestry, rho_estimate, boot_low, boot_high) %>%
  rename(`Spearman's rho` = rho_estimate, `2.5%` = boot_low, `97.5%` = boot_high) %>%
  arrange(desc(feature), desc(group), ancestry)

# print table to file for estimates of spearman's rank correlation
print(xtable(rho,
             digits = c(1, 1, 1, 1, 2, 2, 2),
             type = "latex",
             latex.environments = NULL),
      include.rownames = F,
      file = file_spearmans_rho_ngsadmix)
print(xtable(rho,
             digits = c(1, 1, 1, 1, 2, 2, 2),
             type = "latex",
             latex.environments = NULL),
      include.rownames = F,
      file = file_spearmans_rho_ngsadmix_tex)

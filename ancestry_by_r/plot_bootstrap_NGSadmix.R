#!/usr/bin/env Rscript

# This script plots the the bootstrap results for
# ancestry across recombination quintiles
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(bedr)

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

dir_anc = snakemake@params[["dir_anc"]]
# dir_anc = "ancestry_by_r/results/local_anc_1cM/HILO_MAIZE55/Ne10000_yesBoot"

# png filenames out:
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
png_r5_local_anc = snakemake@output[["png_r5_local_anc"]]
# png_r5_local_anc = "ancestry_by_r/plots/local_anc_by_r_sympatric_only.png"
png_cd5_local_anc = snakemake@output[["png_cd5_local_anc"]]
# png_cd5_local_anc = "ancestry_by_r/plots/local_anc_by_cd_sympatric_only.png"
png_boot_cor_local_anc_bp = snakemake@output[["png_boot_cor_local_anc_bp"]]
# png_boot_cor_local_anc_bp = "ancestry_by_r/plots/boot_cor_local_anc_bp.png"
png_boot_cor_local_anc_r5 = snakemake@output[["png_boot_cor_local_anc_r5"]]
# png_boot_cor_local_anc_r5 = "ancestry_by_r/plots/boot_cor_local_anc_r5.png"


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

# local ancestry by windows:
anc_by_wind = do.call(rbind, lapply(meta_symp$popN, function(i)
  read.table(paste0(dir_anc, "/pop", i, ".anc.wind"),
             stringsAsFactors = F) %>%
    data.table::setnames(c("window", "pop", "anc")))) %>%
  dplyr::left_join(., windows, by = "window") %>%
  dplyr::left_join(., meta_symp, by = "pop") %>%
  dplyr::mutate(., cM_codingMb = 1/(coding_bp/10^5))

# summarise mean across populations:
anc_by_wind_and_zea <- anc_by_wind %>%
  dplyr::mutate(length = end - start) %>%
  dplyr::group_by(zea, window, bin_cd5, bin_r5, length, cM_Mb, cM_codingMb, coding_bp) %>%
  dplyr::summarise(anc = sum(anc*n_local_ancestry)/sum(n_local_ancestry))


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
  guides(color = guide_legend("Subspecies"),
         shape = guide_legend("Subspecies")) +
  ggtitle("Mexicana ancestry by recombination rate") +
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
  guides(color = guide_legend("Subspecies"),
         shape = guide_legend("Subspecies")) +
  ggtitle("Mexicana ancestry by gene density") +
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
  ggtitle("Log(cM/Mb) vs. coding bp in 1cM windows")
ggsave(file = png_cor_r_cd,
       plot = p_cor_r_cd,
       device = "png",
       width = 7, height = 4, 
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
  ggtitle("Mexicana ancestry by recombination rate and elevation") +
  scale_color_manual(values = col_maize_mex_parv) +
  scale_fill_manual(values = col_maize_mex_parv) +
  theme_light() +
  facet_wrap(~bin) +
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
  ggtitle("Mexicana ancestry by recombination rate and elevation") +
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

# mexicana ancestry for K=2 in maize and mexicana, colored by elevation:
p_r5_color_elev <- r5$anc_group_confidence_intervals %>%
  dplyr::filter(ancestry == "mexicana_ancestry") %>%
  dplyr::filter(symp_allo == "sympatric") %>%
  ggplot(aes(x = bin, y = p, group = zea)) +
  # first plot original point estimates for ind. ancestry
  geom_point(data = filter(r5$anc_ind, zea != "parviglumis" & ancestry == "mexicana"),
             aes(x = bin,
                 y = p,
                 shape = zea,
                 color = ELEVATION),
             position = position_jitter(0.2)) +
  scale_color_viridis(direction = -1, option = "viridis") +
  # then add mean for that group
  geom_point(pch = 18, size = .1, alpha = 0.75) +
  # and errorbars for 90% CI around that mean
  # based on bootstrap with NGSAdmix
  geom_errorbar(aes(ymin = low,
                    ymax = high),
                width = .5) +
  xlab("Recombination rate quintile (cM/Mb)") +
  ylab("Proportion mexicana ancestry") +
  theme_classic() +
  facet_grid(zea ~ symp_allo) +
  guides(color = guide_legend("Elevation"),
         shape = guide_legend("Subspecies")) +
  ggtitle("Mexicana ancestry by recombination rate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_r5_color_elev
#ggsave(file = png_r5_color_elev,
#       plot = p_r5_color_elev,
#       device = "png",
#       width = 7,5, height = 4, 
#       units = "in", dpi = 300)



# local ancestry plots:
# recombination rate
# local ancestry plots:
p_r5_local_anc <- anc_by_wind_and_zea %>%
  ggplot(., aes(x = bin_r5, y = anc, group = zea)) +
  # first plot original point estimates for ind. ancestry
  geom_point(aes(shape = zea,
                 color = zea),
             position = position_jitter(0.2)) +
  scale_color_manual(values = col_maize_mex_parv) +
  # then add mean for that group
  geom_point(data = anc_by_wind_and_zea %>%
               dplyr::group_by(bin_r5, zea) %>%
               dplyr::summarise(anc = sum(length*anc)/sum(length)), pch = 18, size = 2) +
  # and errorbars for 90% CI around that mean
  # based on bootstrap
  #geom_errorbar(aes(ymin = low,
  #                  ymax = high),
  #              width = .5) +
  xlab("Recombination rate quintile (cM/Mbp)") +
  ylab("Proportion mexicana ancestry (HMM)") +
  theme_classic() +
  guides(color = guide_legend("Subspecies"),
         shape = guide_legend("Subspecies")) +
  ggtitle("Local mexicana ancestry by recombination rate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p_r5_local_anc
ggsave(file = png_r5_local_anc,
       plot = p_r5_local_anc,
       device = "png",
       width = 5, height = 4, 
       units = "in", dpi = 300)

# coding bp/cM
p_cd5_local_anc <- anc_by_wind_and_zea %>%
  ggplot(., aes(x = bin_cd5, y = anc, group = zea)) +
  # first plot original point estimates for ind. ancestry
  geom_point(aes(shape = zea,
                 color = zea),
             position = position_jitter(0.2)) +
  scale_color_manual(values = col_maize_mex_parv) +
  # then add mean for that group
  geom_point(data = anc_by_wind_and_zea %>%
               dplyr::group_by(bin_cd5, zea) %>%
               dplyr::summarise(anc = sum(length*anc)/sum(length)), pch = 18, size = 2) +
  # and errorbars for 90% CI around that mean
  # based on bootstrap
  #geom_errorbar(aes(ymin = low,
  #                  ymax = high),
  #              width = .5) +
  xlab("Coding density quintile (bp/cM)") +
  ylab("Proportion mexicana ancestry (HMM)") +
  theme_classic() +
  guides(color = guide_legend("Subspecies"),
         shape = guide_legend("Subspecies")) +
  ggtitle("Local mexicana ancestry by coding bp") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p_cd5_local_anc
ggsave(file = png_cd5_local_anc,
       plot = p_cd5_local_anc,
       device = "png",
       width = 5, height = 4, 
       units = "in", dpi = 300)

p_log10_cM_Mb_local_anc <- anc_by_wind_and_zea %>%
  left_join(., windows_by_inv, by = "window") %>%
  left_join(., windows[ , 1:4], by = "window") %>%
  ggplot(., aes(x = log10(cM_Mb), y = anc, group = zea, color = inv)) +
  # first plot original point estimates for ind. ancestry
  geom_point(alpha = 0.5) +
  #scale_color_manual(values = col_maize_mex_parv) +
  xlab("Recombination rate log10(cM/Mb)") +
  ylab("Proportion mexicana ancestry (HMM)") +
  theme_classic() +
  #guides(color = guide_legend("Subspecies"),
  #       shape = guide_legend("Subspecies")) +
  ggtitle("Local mexicana ancestry by recombination rate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~zea)
#p_log10_cM_Mb_local_anc
# NOTE: I SHOULD ADD INV4M TO THIS PLOT

# ------------------------------------------------- #

# boot_local_r5_byBin <- do.call(bind_rows, lapply(0:100, function(i)
#   anc_by_wind_and_zea %>%
#   dplyr::group_by(zea, bin_r5) %>% # re-sample within bins
#   dplyr::sample_n(., size = n(), replace = (i != 0)) %>% # i = 0 is the original sample (sample all w/out replacement)
#     # while i >= 1 is a bootstrap sample with replacement
#   dplyr::group_by(zea) %>%
#   dplyr::summarise(cor.pearson = cor(cM_Mb, anc, method = "pearson"),
#                    cor.log10 = cor(log10(cM_Mb), anc, method = "pearson"),
#                    cor.spearman = cor(cM_Mb, anc, method = "spearman"),
#                    cor.bins_spearman = cor(as.integer(bin_r5), anc, method = "spearman"),
#                    cor.bins_pearson = cor(as.integer(bin_r5), anc, method = "pearson"),
#                    mean_wind = mean(as.integer(bin_r5))) %>%
#       dplyr::mutate(boot = i)))
# # re-sample across all windows, regardless of bin
# test <- data.frame(x = c(1:10),
#                    y = LETTERS[1:10])
# dplyr::sample_frac(test, replace = TRUE)

# plot correlations for local ancestry bootstraps
# color bar for different correlation types
corr_colors = viridis::viridis(5) 
names(corr_colors) = paste0("cor.", c("bins_pearson", "bins_spearman",
                                      "log10", "pearson", "spearman"))
corr_labels_r5 = c("Pearson's corr. anc ~ quintiles", 
                   "Spearman's rank corr. anc ~ quintiles",
                   "Pearson's corr. anc ~ log10(rate)", 
                   "Pearson's corr. anc ~ rate",
                   "Spearman's rank corr. anc ~ rate")
corr_labels_cd5 = c("Pearson's corr. anc ~ quintiles", 
                    "Spearman's rank corr. anc ~ quintiles",
                    "Pearson's corr. anc ~ bp",
                    "Spearman's rank corr. anc ~ bp")

# bootstrap of local ancestry results. resample across all windows,
# regardless of recombination bin
boot_local_r5_all <- do.call(bind_rows, lapply(0:100, function(i)
  anc_by_wind_and_zea %>%
    ungroup(.) %>%
    #left_join(., windows_by_inv, by = "window") %>%
    #filter(inv == "None") %>% # exclude known inversions
    dplyr::sample_n(., size = n(), replace = (i != 0)) %>% # i = 0 is the original sample (sample all w/out replacement)
    # while i >= 1 is a bootstrap sample with replacement
    dplyr::group_by(zea) %>%
    dplyr::summarise(cor.pearson = cor(cM_Mb, anc, method = "pearson"),
                     cor.log10 = cor(log10(cM_Mb), anc, method = "pearson"),
                     cor.spearman = cor(cM_Mb, anc, method = "spearman"),
                     cor.bins_spearman = cor(as.integer(bin_r5), anc, method = "spearman"),
                     cor.bins_pearson = cor(as.integer(bin_r5), anc, method = "pearson"),
                     mean_wind = mean(as.integer(bin_r5))) %>%
    dplyr::mutate(boot = i)))


# plot bootstraps for local ancestry inference
# first get bootstrap confidence intervals
# group means (by symp/allo and zea subspecies, bootstrap = 0 is all the original data)
estimate_r5 <- boot_local_r5_all %>%
  dplyr::filter(boot == 0) %>%
  tidyr::pivot_longer(cols = starts_with("cor."),
                      values_to = "correlation",
                      names_to = "method")

# get percentiles from bootstrap
boot_perc_r5 <- boot_local_r5_all %>%
  dplyr::filter(boot != 0) %>% # exclude boot = 0, which is the original sample
  tidyr::pivot_longer(cols = starts_with("cor."),
                      values_to = "correlation",
                      names_to = "method") %>%
  dplyr::group_by(method, zea) %>%
  dplyr::summarise(low_boot = quantile(correlation, alpha/2), # low and high bounds
                   high_boot = quantile(correlation, 1-alpha/2), # of 90% conf. interval
                   median_boot = median(correlation),
                   mean_boot = mean(correlation))

# calculate basic bootstrap (i.e. 'pivot') confidence intervals together (to plot later as range)
boot_confidence_intervals_r5 <- estimate_r5 %>%
  left_join(., boot_perc_r5, by = c("zea", "method")) %>%
  mutate(low = 2*correlation - high_boot, # p is the estimate from the original sample
         high = 2*correlation - low_boot)

p_boot_cor_local_anc_r5 <- boot_local_r5_all %>%
  tidyr::pivot_longer(cols = starts_with("cor."),
                      values_to = "correlation",
                      names_to = "method") %>%
  dplyr::filter(boot != 0) %>% # plot only bootstraps, not original
  #mutate(method = sapply(.$method, function(x) substr(method, 5, Inf))) %>%
  ggplot(., aes(x = method, y = correlation, fill = method, color = method)) +
  theme_light() +
  geom_violin(alpha = 0.5, lwd = 0.2) + 
  facet_wrap(~zea) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Proportion mexicana local ancestry by recombination rate (95% CI simple bootstrap)") +
  scale_fill_manual(values = corr_colors, labels = corr_labels_r5) +
  scale_color_manual(values = corr_colors, labels = corr_labels_r5) +
  geom_errorbar(data = boot_confidence_intervals_r5,
                aes(ymin = low, ymax = high), color = "black") +
  geom_point(data = boot_confidence_intervals_r5, color = "black")
p_boot_cor_local_anc_r5
ggsave(png_boot_cor_local_anc_r5,
       plot = p_boot_cor_local_anc_r5,
       device = "png", dpi = 300,
       height = 5, width = 8, units = "in")

# plot bootstrapped correlation estimates for cM/(coding bp):
boot_local_cd5_all <- do.call(bind_rows, lapply(0:100, function(i)
  anc_by_wind_and_zea %>% # note, can't use cM/(coding_Mb) because for some windows that's a dividing by zero error
    ungroup(.) %>%
    dplyr::sample_n(., size = n(), replace = (i != 0)) %>% # i = 0 is the original sample (sample all w/out replacement)
    # while i >= 1 is a bootstrap sample with replacement
    dplyr::group_by(zea) %>%
    dplyr::summarise(cor.pearson = cor(coding_bp, anc, method = "pearson"),
                     #cor.log10 = cor(log10(coding_bp), anc, method = "pearson"), # no log10 b/c log10(0) is -Inf
                     cor.spearman = cor(coding_bp, anc, method = "spearman"),
                     cor.bins_spearman = cor(as.integer(bin_cd5), anc, method = "spearman"),
                     cor.bins_pearson = cor(as.integer(bin_cd5), anc, method = "pearson"),
                     mean_wind = mean(as.integer(bin_cd5))) %>%
    dplyr::mutate(boot = i)))


# plot bootstraps for local ancestry inference
# first get bootstrap confidence intervals
# group means (by symp/allo and zea subspecies, bootstrap = 0 is all the original data)
estimate_cd5 <- boot_local_cd5_all %>%
  dplyr::filter(boot == 0) %>%
  tidyr::pivot_longer(cols = starts_with("cor."),
                      values_to = "correlation",
                      names_to = "method")

# get percentiles from bootstrap
boot_perc_cd5 <- boot_local_cd5_all %>%
  dplyr::filter(boot != 0) %>% # exclude boot = 0, which is the original sample
  tidyr::pivot_longer(cols = starts_with("cor."),
                      values_to = "correlation",
                      names_to = "method") %>%
  dplyr::group_by(method, zea) %>%
  dplyr::summarise(low_boot = quantile(correlation, alpha/2), # low and high bounds
                   high_boot = quantile(correlation, 1-alpha/2), # of 90% conf. interval
                   median_boot = median(correlation),
                   mean_boot = mean(correlation))

# calculate basic bootstrap (i.e. 'pivot') confidence intervals together (to plot later as range)
boot_confidence_intervals_cd5 <- estimate_cd5 %>%
  left_join(., boot_perc_cd5, by = c("zea", "method")) %>%
  mutate(low = 2*correlation - high_boot, # p is the estimate from the original sample
         high = 2*correlation - low_boot)

p_boot_cor_local_anc_cd5 <- boot_local_cd5_all %>%
  tidyr::pivot_longer(cols = starts_with("cor."),
                      values_to = "correlation",
                      names_to = "method") %>%
  dplyr::filter(boot != 0) %>% # plot only bootstraps, not original
  #mutate(method = sapply(.$method, function(x) substr(method, 5, Inf))) %>%
  ggplot(., aes(x = method, y = correlation, fill = method, color = method)) +
  theme_light() +
  geom_violin(alpha = 0.5, lwd = 0.2) + 
  facet_wrap(~zea) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Proportion mexicana local ancestry by coding bp in 1cM window (95% CI simple bootstrap)") +
  scale_fill_manual(values = corr_colors, labels = corr_labels_cd5) +
  scale_color_manual(values = corr_colors, labels = corr_labels_cd5) +
  geom_errorbar(data = boot_confidence_intervals_cd5,
                aes(ymin = low, ymax = high), color = "black") +
  geom_point(data = boot_confidence_intervals_cd5, color = "black")
p_boot_cor_local_anc_cd5

# plot ancestry correlations with all bp and coding bp (at 1cM resolution) together
boot_local_all <- do.call(bind_rows, lapply(0:100, function(i)
  anc_by_wind_and_zea %>%
    #left_join(., windows_by_inv, by = "window") %>%
    #filter(inv == "None") %>% # exclude known inversions
    ungroup(.) %>%
    dplyr::sample_n(., size = n(), replace = (i != 0)) %>% # i = 0 is the original sample (sample all w/out replacement)
    # while i >= 1 is a bootstrap sample with replacement
    mutate(bp_r5 = 1/(cM_Mb)*10^6) %>% # all bp in 1cM window
    rename(bp_cd5 = coding_bp) %>%
    tidyr::pivot_longer(data = ., cols = starts_with("bp"),
                        names_to = "bp_type", values_to = "bp") %>%
    mutate(count_type = ifelse(bp_type == "bp_r5", "all", "coding"),
           bin = ifelse(count_type == "all", 
                        6 - as.integer(bin_r5), # reverse order of bins for r5 b/c higher recombination rates = low # bp in a 1cM window
                        as.integer(bin_cd5))) %>%
    dplyr::group_by(zea, count_type) %>%
    dplyr::summarise(cor.pearson = cor(bp, anc, method = "pearson"),
                     #cor.log10 = cor(log10(coding_bp), anc, method = "pearson"), # no log10 b/c log10(0) is -Inf
                     cor.spearman = cor(bp, anc, method = "spearman"),
                     cor.bins_spearman = cor(bin, anc, method = "spearman"),
                     cor.bins_pearson = cor(bin, anc, method = "pearson"),
                     mean_wind = mean(bin)) %>%
    dplyr::mutate(boot = i)))


# plot bootstraps for local ancestry inference
# first get bootstrap confidence intervals
# group means (by symp/allo and zea subspecies, bootstrap = 0 is all the original data)
estimate <- boot_local_all %>%
  dplyr::filter(boot == 0) %>%
  tidyr::pivot_longer(cols = starts_with("cor."),
                      values_to = "correlation",
                      names_to = "method")

# get percentiles from bootstrap
boot_perc <- boot_local_all %>%
  dplyr::filter(boot != 0) %>% # exclude boot = 0, which is the original sample
  tidyr::pivot_longer(cols = starts_with("cor."),
                      values_to = "correlation",
                      names_to = "method") %>%
  dplyr::group_by(method, count_type, zea) %>%
  dplyr::summarise(low_boot = quantile(correlation, alpha/2), # low and high bounds
                   high_boot = quantile(correlation, 1-alpha/2), # of 90% conf. interval
                   median_boot = median(correlation),
                   mean_boot = mean(correlation))

# calculate basic bootstrap (i.e. 'pivot') confidence intervals together (to plot later as range)
boot_confidence_intervals <- estimate %>%
  left_join(., boot_perc, by = c("zea", "method", "count_type")) %>%
  mutate(low = 2*correlation - high_boot, # p is the estimate from the original sample
         high = 2*correlation - low_boot)

p_boot_cor_local_anc_bp <- boot_local_all %>%
  tidyr::pivot_longer(cols = starts_with("cor."),
                      values_to = "correlation",
                      names_to = "method") %>%
  dplyr::filter(boot != 0) %>% # plot only bootstraps, not original
  ggplot(., aes(x = method, y = correlation, fill = method, color = method)) +
  theme_light() +
  geom_violin(alpha = 0.5, lwd = 0.2) + 
  facet_grid(zea ~ count_type) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Proportion mexicana local ancestry by bp in 1cM window (95% CI simple bootstrap)") +
  scale_fill_manual(values = corr_colors, labels = corr_labels_cd5) +
  scale_color_manual(values = corr_colors, labels = corr_labels_cd5) +
  geom_errorbar(data = boot_confidence_intervals,
                aes(ymin = low, ymax = high), color = "black") +
  geom_point(data = boot_confidence_intervals, color = "black")
p_boot_cor_local_anc_bp
ggsave(png_boot_cor_local_anc_bp,
       plot = p_boot_cor_local_anc_bp,
       device = "png", dpi = 300,
       height = 5, width = 8, units = "in")


# ---------------------------------------------------- #
# # get data for supplement K=3 figures
# anc_boot_3 <- do.call(rbind,
#                lapply(0:100, function(BOOT) do.call(rbind,
#                                           lapply(1:5, function(r)
#                                             read.table(paste0("results/bootstrap/windows_", WIND,
#                          "/r5_recomb", r, "/", PREFIX_3, "/K", K_3, "/boot", BOOT, ".anc"),
#                   header = T, sep = "\t") %>%
#   mutate(bootstrap = BOOT) %>%
#   mutate(r5_bin = r))))) %>%
#   left_join(., meta, by = "ID")
# 
# # set sequencing coverage cutoff for inclusion
# min_coverage = 0.1
# 
# anc_boot_mean_3 <- anc_boot_3 %>%
#   filter(zea == "parviglumis" |
#            est_coverage >= min_coverage) %>%
#   group_by(group, zea, symp_allo, bootstrap, r5_bin) %>%
#   summarise(mexicana_ancestry = mean(mexicana),
#             maize_ancestry = mean(maize),
#             parviglumis_ancestry = mean(parviglumis))
# # point estimate individual ancestry
# anc_ind_mean_3 <- anc_boot_3 %>%
#   filter(zea == "parviglumis" |
#            est_coverage >= min_coverage) %>%
#   filter(bootstrap == 1) %>%
#   gather(., key = "ancestry", value = "p", ancestries)
# # bootstrap 90% confidence interval individual ancestry
# anc_ind_boot_confidence_intervals_3 <- anc_boot_3 %>%
#   filter(bootstrap != 0) %>% # this is the original sample
#   gather(., key = "ancestry", value = "p",
#          ancestries) %>%
#   group_by(ID, ancestry, r5_bin) %>%
#   summarise(low_boot = nth(p, 6, order_by = p), # low and high bounds
#             high_boot = nth(p, 95, order_by = p), # of 90% conf. interval
#             mid_boot = nth(p, 50, order_by = p),
#             mean_boot = mean(p))
# anc_ind_estimate_3 <- anc_ind_mean_3 %>%
#   left_join(., anc_ind_boot_confidence_intervals_3, by = c("ancestry", "ID", "r5_bin")) %>%
#   mutate(low = low_boot - mean_boot + p,
#          high = high_boot - mean_boot + p)
# 
# 
# # plot
# anc_boot_mean_3 %>%
#   ggplot(., aes(group = r5_bin, y = maize_ancestry)) +
#   geom_boxplot() +
#   facet_wrap(~group) +
#   ggtitle("Mean maize ancestry K3 bootstrapped")
# anc_boot_mean_3 %>%
#   ggplot(., aes(group = r5_bin, y = mexicana_ancestry)) +
#   geom_boxplot() +
#   facet_wrap(~group) +
#   ggtitle("Mean mexicana ancestry K3 bootstrapped")
# anc_boot_mean_3 %>%
#   ggplot(., aes(group = r5_bin, y = parviglumis_ancestry)) +
#   geom_boxplot() +
#   facet_wrap(~group) +
#   ggtitle("Mean parviglumis ancestry K3 bootstrapped")
# 
# # plot all together
# anc_boot_mean_3 %>%
#   gather(., key = "ancestry", value = "p",
#          paste(ancestries, "ancestry", sep = "_")) %>%
#   ggplot(aes(x = as.factor(r5_bin), y = p, color = ancestry)) +
#   geom_boxplot() +
#   facet_wrap(~group) +
#   xlab("recombination quintile (low to high r)") +
#   ggtitle("Mean K3 ancestries bootstrapped")
# 
# anc_boot_mean_3 %>%
#   ggplot(aes(x = as.factor(r5_bin), y = mexicana_ancestry,
#              color = group)) +
#   geom_boxplot(outlier.size = 1, notch = T) +
#   geom_jitter(shape = 16, size = .1, position=position_jitter(0.2)) +
#   #geom_point(size = .1, color = "black") +
#   facet_wrap(~group) +
#   ggtitle("Mean K3 mexicana ancestry bootstrapped")
# 
# anc_group_mean_3 <- anc_boot_mean_3 %>%
#   filter(bootstrap == 0) %>%
#   gather(., key = "ancestry", value = "p",
#          paste(ancestries, "ancestry", sep = "_"))
# anc_boot_confidence_intervals_3 <- anc_boot_mean_3 %>%
#   filter(bootstrap != 0) %>% # this is the original sample
#   gather(., key = "ancestry", value = "p",
#          paste(ancestries, "ancestry", sep = "_")) %>%
#   group_by(ancestry, group, r5_bin) %>%
#   summarise(low_boot = nth(p, 6, order_by = p), # low and high bounds
#             high_boot = nth(p, 95, order_by = p), # of 90% conf. interval
#             mid_boot = nth(p, 50, order_by = p),
#             mean_boot = mean(p))
# anc_group_estimate_3 <- anc_group_mean_3 %>%
#   left_join(., anc_boot_confidence_intervals_3, by = c("ancestry", "group", "r5_bin")) %>%
#   mutate(low = low_boot - mean_boot + p,
#          high = high_boot - mean_boot + p)
# 
# # plot 90% confidence interval around mean from bootstrap:
# anc_group_estimate_3 %>%
#   ggplot(aes(x = as.factor(r5_bin), y = p,
#              color = ancestry)) +
#   geom_errorbar(aes(ymin = low,
#                     ymax = high),
#                 width = .5) +
#   geom_point(pch = 18, size = 1) +
#   facet_wrap(~group) +
#   scale_color_manual(values = col_maize_mex_parv) +
#   #theme_classic() +
#   xlab("recombination rate quintile (low to high r)") +
#   ggtitle("Mean K=3 ancestry estimates and 90% bootstrap CI")
# ggsave("plots/K3_mean_ancestry_all_bootstrap_90CI.png",
#        width = 8, height = 8, units = "in")
# 
# # look at just mexicana ancestry in maize and mexicana w/in K=3
# anc_group_estimate_3 %>%
#   filter(zea != "parviglumis") %>%
#   filter(ancestry == "mexicana_ancestry") %>%
#   ggplot(aes(x = as.factor(r5_bin), y = p, group = zea)) +
#   # first plot original point estimates for ind. ancestry
#   geom_point(data = filter(anc_ind_mean_3, zea != "parviglumis" & ancestry == "mexicana") %>%
#                mutate(depth = ifelse(est_coverage > 2, 2, est_coverage)),
#              aes(x = as.factor(r5_bin),
#                  y = p,
#                  #size = depth
#                  color = zea),
#              shape = 16,
#              position = position_jitter(0.2)) +
#   scale_color_manual(values = col_maize_mex_parv) +
#   #scale_size(range = c(.05, 2)) +
#   # then add mean for that group
#   geom_point(pch = 18, size = 2) +
#   # and errorbars for 90% CI around that mean
#   # based on bootstrap with NGSAdmix
#   geom_errorbar(aes(ymin = low,
#                     ymax = high),
#                 width = .5) +
#   facet_wrap(~symp_allo) +
#   xlab("recombination rate quintile (low to high r)") +
#   ylab("proportion mexicana ancestry") +
#   theme(legend.title = element_blank()) +
#   ggtitle("Mexicana ancestry by recombination rate")
# ggsave("plots/K3_ind_ancestries_jitter_and_mean_bootstrap_90CI.png",
#        width = 8, height = 8, units = "in")
# 
# # what are the spearman's rank correlations for NGSAdmix mexicana ancestry estimates vs. recombination rate quintile?
# filter(r5$anc_boot_mean, zea != "parviglumis" & symp_allo == "sympatric") %>%
#   group_by(zea, bootstrap) %>%
#   filter(bootstrap != 0) %>%
#   summarise(spearmans.cor = cor(quintile, mexicana_ancestry, method = "spearman")) %>%
#   group_by(zea) %>%
#   summarise(#min = min(spearmans.cor),
#             median = median(spearmans.cor),
#             #max = max(spearmans.cor),
#             q_.025 = quantile(spearmans.cor, 0.025),
#             q_.975 = quantile(spearmans.cor, 0.975)) %>%
#   left_join(filter(r5$anc_boot_mean, zea != "parviglumis" & symp_allo == "sympatric") %>%
#               filter(bootstrap == 0) %>%
#               group_by(zea) %>%
#               summarise(rho = cor(quintile, mexicana_ancestry, method = "spearman")),
#             .,
#             by = "zea") %>%
#   mutate(., lower = 2*rho - q_.975,
#          upper = 2*rho - q_.025)
# # vs. coding density quintile?
# filter(cd5$anc_boot_mean, zea != "parviglumis" & symp_allo == "sympatric") %>%
#   group_by(zea, bootstrap) %>%
#   filter(bootstrap != 0) %>%
#   summarise(spearmans.cor = cor(quintile, mexicana_ancestry, method = "spearman")) %>%
#   group_by(zea) %>%
#   summarise(median = median(spearmans.cor),
#     q_.025 = quantile(spearmans.cor, 0.025),
#     q_.975 = quantile(spearmans.cor, 0.975)) %>%
#   left_join(filter(cd5$anc_boot_mean, zea != "parviglumis" & symp_allo == "sympatric") %>%
#               filter(bootstrap == 0) %>%
#               group_by(zea) %>%
#               summarise(rho = cor(quintile, mexicana_ancestry, method = "spearman")),
#             .,
#             by = "zea") %>%
#   mutate(., lower = 2*rho - q_.975,
#          upper = 2*rho - q_.025)

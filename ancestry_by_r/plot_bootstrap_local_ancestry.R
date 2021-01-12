#!/usr/bin/env Rscript

# this script makes local ancestry plots by 
# recombination rate and gene density
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

source(snakemake@input[["colors"]]) # plotting colors
# source("colors.R")
load(snakemake@input[["meta"]]) # meta
# load("samples/HILO_MAIZE55_meta.RData")

dir_anc = snakemake@params[["dir_anc"]]
# dir_anc = "ancestry_by_r/results/local_anc_1cM/HILO_MAIZE55/Ne10000_yesBoot"

# png filenames out:
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





# recombination rate local ancestry plots:
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



rho = bind_rows(mutate(r5$spearman, 
                       method = "NGSAdmix",
                       feature = "recombination rate (cM/Mb)"),
                mutate(cd5$spearman,
                       method = "NGSAdmix",
                       feature = "gene density (coding bp/cM)")) %>%
  mutate(group = stringr::str_replace(group, "_", " "),
         resolution = "genomic quintiles") %>%
  dplyr::select(method, feature, resolution, group, rho_estimate, boot_low, boot_high) %>%
  rename(`Pearson's rank correlation` = rho_estimate, `2.5%` = boot_low, `97.5%` = boot_high)

# print table to file for estimates of Pearson's rank correlation
print(xtable(rho, 
             digits = 3,
             label = "tbl_pearsons_rho_ngsadmix",
             type = "latex", 
             latex.environments = NULL),
      include.rownames = F,
      file = file_pearsons_rho_ngsadmix)



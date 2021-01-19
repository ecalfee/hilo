#!/usr/bin/env Rscript

# this script makes local ancestry plots by 
# recombination rate and gene density
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(bedr)
library(xtable)
library(boot)

# load variables from Snakefile
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
n_boot = snakemake@params[["n_boot"]] # number of bootstraps
# n_boot = 10000
alpha = 0.05 # use 95% confidence level for bootstraps

# png and table filenames out:
png_r5_local_anc = snakemake@output[["png_r5_local_anc"]]
# png_r5_local_anc = "ancestry_by_r/plots/local_anc_by_r_quintiles.png"
png_r_cont_local_anc = snakemake@output[["png_r_cont_local_anc"]]
# png_r_cont_local_anc = "ancestry_by_r/plots/local_anc_by_r_continuous.png"
png_cd5_local_anc = snakemake@output[["png_cd5_local_anc"]]
# png_cd5_local_anc = "ancestry_by_r/plots/local_anc_by_cd_quintiles.png"
file_pearsons_rho = snakemake@output[["file_pearsons_rho_local_anc"]]
# file_pearsons_rho = "ancestry_by_r/tables/pearsons_rho_local_ancestry.tex"

random_seed = 7135891 # for reproducing bootstrap results

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
  dplyr::group_by(zea, window, bin_cd5, bin_r5, quintile_cd5, quintile_r5, length, cM_Mb, cM_codingMb, coding_bp) %>%
  dplyr::summarise(anc = sum(anc*n_local_ancestry)/sum(n_local_ancestry),
                   .groups = "drop")

# bootstrap across windows:
# re-sampling with replacement is done across all windows
# and calculates spearman's rank correlations as well as mean ancestry per quintile
# for both recombination rate and coding density
calc_spearmans_cor <- function(d, i) { # d is data (only 1 zea), i is indices
  
  spearman = d[i, ] %>%
    dplyr::summarise(cor.spearman.r = cor(cM_Mb, anc, method = "spearman"),
                     cor.spearman.cd = cor(coding_bp, anc, method = "spearman"))

  r_quintiles = d[i, ] %>%
    dplyr::group_by(bin_r5, quintile_r5) %>%
    dplyr::summarise(anc = sum(length*anc)/sum(length), # take weighted mean by window size (bp), so comparable to NGSAdmix and f4 estimates by quintile
                     .groups = "drop") %>%
    pivot_wider(., 
                names_from = quintile_r5, 
                names_prefix = "r", 
                values_from = anc, 
                id_cols = quintile_r5)
  
  cd_quintiles = d[i, ] %>%
    dplyr::group_by(bin_cd5, quintile_cd5) %>%
    dplyr::summarise(anc = sum(length*anc)/sum(length), .groups = "drop") %>%
    pivot_wider(., 
                names_from = quintile_cd5, 
                names_prefix = "cd", 
                values_from = anc, 
                id_cols = quintile_cd5)
  
  return(unlist(bind_cols(spearman, r_quintiles, cd_quintiles)))
}

# bootstrap resampling each recombination rate quintile independently
set.seed(random_seed)
boot_maize <- boot(data = filter(anc_by_wind_and_zea, zea == "maize"), 
                statistic = calc_spearmans_cor,
                sim = "ordinary",
                stype = "i",
                R = n_boot)
set.seed(random_seed)
boot_mexicana <- boot(data = filter(anc_by_wind_and_zea, zea == "mexicana"), 
                   statistic = calc_spearmans_cor,
                   sim = "ordinary",
                   stype = "i",
                   R = n_boot)

# get 95% percentile intervals
ci_boot_maize <- t(sapply(1:length(boot_maize$t0), function(x)
  boot.ci(boot_maize, index = x, conf = 0.95, type = "perc")$perc[4:5])) %>%
  data.frame(.) %>%
  data.table::setnames(c("low", "high")) %>%
  mutate(feature = names(boot_maize$t0),
         estimate = boot_maize$t0,
         zea = "maize")
ci_boot_mexicana <- t(sapply(1:length(boot_mexicana$t0), function(x)
  boot.ci(boot_mexicana, index = x, conf = 0.95, type = "perc")$perc[4:5])) %>%
  data.frame(.) %>%
  data.table::setnames(c("low", "high")) %>%
  mutate(feature = names(boot_mexicana$t0),
         estimate = boot_mexicana$t0,
         zea = "mexicana")

# all confidence intervals
ci_boot <- bind_rows(ci_boot_maize, ci_boot_mexicana)

# confidence intervals for recombination quintiles in tidy form for plotting
ci_r5 <- ci_boot %>%
  filter(feature %in% paste0("r", 1:5)) %>%
  dplyr::mutate(quintile_r5 = as.integer(substr(feature, 2, 3))) %>%
  left_join(., anc_by_wind_and_zea %>%
              dplyr::group_by(bin_r5, quintile_r5, zea) %>%
              dplyr::summarise(anc = sum(length*anc)/sum(length), .groups = "drop"),
            by = c("zea", "quintile_r5"))

# confidence intervals for coding density quintiles in tidy form for plotting
ci_cd5 <- ci_boot %>%
  filter(feature %in% paste0("cd", 1:5)) %>%
  dplyr::mutate(quintile_cd5 = as.integer(substr(feature, 3, 4))) %>%
  left_join(., anc_by_wind_and_zea %>%
              dplyr::group_by(bin_cd5, quintile_cd5, zea) %>%
              dplyr::summarise(anc = sum(length*anc)/sum(length), .groups = "drop"),
            by = c("zea", "quintile_cd5"))

# recombination rate local ancestry plots:
p_r5_local_anc <- anc_by_wind_and_zea %>%
  ggplot(., aes(x = bin_r5, y = anc, group = zea)) +
  # first plot original point estimates for ind. ancestry
  geom_point(aes(shape = zea,
                 color = zea),
             alpha = 0.75,
             position = position_jitter(0.2)) +
  scale_color_manual(values = col_maize_mex_parv) +
  # then add mean for that group
  geom_point(data = ci_r5,
             aes(y = estimate),
             pch = 18, size = 2) +
  # and errorbars for 95% CI around that mean
  # based on bootstrap
  geom_errorbar(
    data = ci_r5,
    aes(ymin = low,
        ymax = high),
    width = .5) +
  xlab("Recombination rate quintile (cM/Mbp)") +
  ylab("Proportion mexicana ancestry (HMM)") +
  theme_classic() +
  #guides(color = guide_legend("Subspecies"),
  #       shape = guide_legend("Subspecies")) +
  #ggtitle("Local mexicana ancestry by recombination rate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None") +
  facet_wrap(~zea)
#p_r5_local_anc
ggsave(file = png_r5_local_anc,
       plot = p_r5_local_anc,
       device = "png",
       width = 7.5, height = 4, 
       units = "in", dpi = 300)

p_r_cont_local_anc <- anc_by_wind_and_zea %>%
  ggplot(., aes(x = log10(cM_Mb), y = anc, group = zea)) +
  # first plot original point estimates for ind. ancestry
  geom_point(aes(shape = zea,
                 color = zea),
             alpha = 0.75) +
  scale_color_manual(values = col_maize_mex_parv) +
  xlab("Recombination rate log10(cM/Mbp)") +
  ylab("Proportion mexicana ancestry (HMM)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None") +
  facet_wrap(~zea) + 
  geom_smooth(method = "lm", color = "black", se = F)
#p_r_cont_local_anc
ggsave(file = png_r_cont_local_anc,
       plot = p_r_cont_local_anc,
       device = "png",
       width = 7.5, height = 4, 
       units = "in", dpi = 300)


# coding bp/cM
p_cd5_local_anc <- anc_by_wind_and_zea %>%
  ggplot(., aes(x = bin_r5, y = anc, group = zea)) +
  # first plot original point estimates for ind. ancestry
  geom_point(aes(shape = zea,
                 color = zea),
             alpha = 0.75,
             position = position_jitter(0.2)) +
  scale_color_manual(values = col_maize_mex_parv) +
  # then add mean for that group
  geom_point(data = ci_r5,
             aes(y = estimate),
             pch = 18, size = 2) +
  # and errorbars for 95% CI around that mean
  # based on bootstrap
  geom_errorbar(
    data = ci_r5,
    aes(ymin = low,
        ymax = high),
    width = .5) +
  xlab("Coding density quintile (bp/cM)") +
  ylab("Proportion mexicana ancestry (HMM)") +
  theme_classic() +
  #guides(color = guide_legend("Subspecies"),
  #       shape = guide_legend("Subspecies")) +
  #ggtitle("Local mexicana ancestry by coding bp") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None") +
  facet_wrap(~zea)
#p_cd5_local_anc
ggsave(file = png_cd5_local_anc,
       plot = p_cd5_local_anc,
       device = "png",
       width = 7.5, height = 4, 
       units = "in", dpi = 300)

# make tidy table of spearman's rank correlation (and bootstrap ci) results
rho <- ci_boot %>%
  filter(feature %in% c("cor.spearman.r", "cor.spearman.cd")) %>%
  mutate(feature = ifelse(feature == "cor.spearman.r", "recombination rate (cM/Mb)",
                          "gene density (coding bp/cM)"),
         method = "ancestry_hmm",
         group = paste("sympatric", zea, sep = "_"),
         resolution = "1cM windows") %>%
  dplyr::select(method, group, feature, resolution, estimate, low, high) %>%
  rename(`Pearson's rank correlation` = estimate, `2.5%` = low, `97.5%` = high)

# print table to file for estimates of Pearson's rank correlation
print(xtable(rho, 
             digits = 4,
             label = "tbl_pearsons_rho_local_ancestry",
             type = "latex", 
             latex.environments = NULL),
      include.rownames = F,
      file = file_pearsons_rho)



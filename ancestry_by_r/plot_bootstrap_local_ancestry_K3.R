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
# load("samples/HILO_MAIZE55_PARV50_meta.RData")
dir_anc = snakemake@params[["dir_anc"]]
# dir_anc = "ancestry_by_r/results/local_anc_1cM/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot"
n_boot = snakemake@params[["n_boot"]] # number of bootstraps
# n_boot = 10000
alpha = 0.05 # use 95% confidence level for bootstraps
ancestries = c("maize", "mexicana", "parv")

# png and table filenames out:
png_r5_local_anc = snakemake@output[["png_r5_local_anc"]]
# png_r5_local_anc = "ancestry_by_r/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_local_anc_by_r_quintiles.png"

png_r_cont_local_anc = snakemake@output[["png_r_cont_local_anc"]]
# png_r_cont_local_anc = "ancestry_by_r/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_local_anc_by_r_continuous.png"
png_r_cont_local_anc_lzw = snakemake@output[["png_r_cont_local_anc_lzw"]]
# png_r_cont_local_anc_lzw = "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_local_anc_by_r_continuous.tif"

png_cd5_local_anc = snakemake@output[["png_cd5_local_anc"]]
# png_cd5_local_anc = "ancestry_by_r/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_local_anc_by_cd_quintiles.png"

file_spearmans_rho = snakemake@output[["file_spearmans_rho_local_anc"]]
# file_spearmans_rho = "ancestry_by_r/tables/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_spearmans_rho_local_ancestry.tex"
file_spearmans_rho_tex = snakemake@output[["file_spearmans_rho_local_anc_tex"]]
# file_spearmans_rho_tex = "../hilo_manuscript/tables/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_spearmans_rho_local_ancestry.tex"

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
  do.call(rbind, lapply(ancestries, function(a)
    read.table(paste0(dir_anc, "/", a, "_anc/pop", i, ".anc.wind"),
               stringsAsFactors = F) %>%
    data.table::setnames(c("window", "pop", "p")) %>%
      dplyr::mutate(ancestry = a))))) %>%
  dplyr::left_join(., windows, by = "window") %>%
  dplyr::left_join(., meta_symp, by = "pop") %>%
  dplyr::mutate(., cM_codingMb = 1/(coding_bp/10^5))

# summarise mean across populations:
anc_by_wind_and_zea <- anc_by_wind %>%
  dplyr::mutate(length = end - start) %>%
  dplyr::group_by(zea, ancestry, window, bin_cd5, bin_r5, quintile_cd5, quintile_r5, length, cM_Mb, cM_codingMb, coding_bp) %>%
  dplyr::summarise(p = sum(p*n_local_ancestry)/sum(n_local_ancestry),
                   .groups = "drop")

# bootstrap across windows:
# re-sampling with replacement is done across all windows
# and calculates spearman's rank correlations as well as mean ancestry per quintile
# for both recombination rate and coding density
calc_spearmans_cor <- function(d, i) { # d is data (only 1 zea), i is indices
  
  spearman = d[i, ] %>%
    dplyr::summarise(cor.spearman.r = cor(cM_Mb, p, method = "spearman"),
                     cor.spearman.cd = cor(coding_bp, p, method = "spearman"))

  r_quintiles = d[i, ] %>%
    dplyr::group_by(bin_r5, quintile_r5) %>%
    dplyr::summarise(p = sum(length*p)/sum(length), # take weighted mean by window size (bp), so comparable to NGSAdmix and f4 estimates by quintile
                     .groups = "drop") %>%
    pivot_wider(., 
                names_from = quintile_r5, 
                names_prefix = "r", 
                values_from = p, 
                id_cols = quintile_r5)
  
  cd_quintiles = d[i, ] %>%
    dplyr::group_by(bin_cd5, quintile_cd5) %>%
    dplyr::summarise(p = sum(length*p)/sum(length), .groups = "drop") %>%
    pivot_wider(., 
                names_from = quintile_cd5, 
                names_prefix = "cd", 
                values_from = p, 
                id_cols = quintile_cd5)
  
  return(unlist(bind_cols(spearman, r_quintiles, cd_quintiles)))
}

# bootstrap resampling each recombination rate quintile independently
# make lists for each ancestry...
boot_maize = list()
boot_mexicana = list()

set.seed(random_seed)
for (a in ancestries){
  boot_maize[[a]] <- boot(data = filter(anc_by_wind_and_zea, zea == "maize", ancestry == a), 
                          statistic = calc_spearmans_cor,
                          sim = "ordinary",
                          stype = "i",
                          R = n_boot)
}

set.seed(random_seed)
for (a in ancestries){
  boot_mexicana[[a]] <- boot(data = filter(anc_by_wind_and_zea, zea == "mexicana", ancestry == a), 
                   statistic = calc_spearmans_cor,
                   sim = "ordinary",
                   stype = "i",
                   R = n_boot)
}
# get 95% percentile intervals
ci_boot_maize = list()
ci_boot_mexicana = list()
for (a in ancestries){
  ci_boot_maize[[a]] <- t(sapply(1:length(boot_maize[[a]]$t0), function(x)
    boot.ci(boot_maize[[a]], index = x, conf = 0.95, type = "perc")$perc[4:5])) %>%
    data.frame(.) %>%
    data.table::setnames(c("low", "high")) %>%
    mutate(feature = names(boot_maize[[a]]$t0),
           estimate = boot_maize[[a]]$t0,
           zea = "maize",
           ancestry = a)
  ci_boot_mexicana[[a]] <- t(sapply(1:length(boot_mexicana[[a]]$t0), function(x)
    boot.ci(boot_mexicana[[a]], index = x, conf = 0.95, type = "perc")$perc[4:5])) %>%
    data.frame(.) %>%
    data.table::setnames(c("low", "high")) %>%
    mutate(feature = names(boot_mexicana[[a]]$t0),
           estimate = boot_mexicana[[a]]$t0,
           zea = "mexicana",
           ancestry = a)
}


# all confidence intervals
ci_boot <- bind_rows(do.call(bind_rows, ci_boot_maize), 
                     do.call(bind_rows, ci_boot_mexicana))

# confidence intervals for recombination quintiles in tidy form for plotting
ci_r5 <- ci_boot %>%
  filter(feature %in% paste0("r", 1:5)) %>%
  dplyr::mutate(quintile_r5 = as.integer(substr(feature, 2, 3))) %>%
  left_join(., anc_by_wind_and_zea %>%
              dplyr::group_by(bin_r5, quintile_r5, zea, ancestry) %>%
              dplyr::summarise(p = sum(length*p)/sum(length), .groups = "drop"),
            by = c("zea", "ancestry", "quintile_r5"))

# confidence intervals for coding density quintiles in tidy form for plotting
ci_cd5 <- ci_boot %>%
  filter(feature %in% paste0("cd", 1:5)) %>%
  dplyr::mutate(quintile_cd5 = as.integer(substr(feature, 3, 4))) %>%
  left_join(., anc_by_wind_and_zea %>%
              dplyr::group_by(bin_cd5, quintile_cd5, zea, ancestry) %>%
              dplyr::summarise(p = sum(length*p)/sum(length), .groups = "drop"),
            by = c("zea", "ancestry", "quintile_cd5"))

# recombination rate local ancestry plots:
p_r5_local_anc <- anc_by_wind_and_zea %>%
  dplyr::mutate(ancestry = ifelse(ancestry == "parv", "parviglumis", ancestry)) %>%
  ggplot(., aes(x = bin_r5, y = p, group = paste(zea, ancestry))) +
  # first plot original point estimates for ind. ancestry
  geom_point(aes(shape = zea,
                 color = ancestry),
             alpha = 0.75,
             position = position_jitter(0.2)) +
  scale_color_manual(values = col_maize_mex_parv) +
  # then add mean for that group
  geom_point(data = ci_r5 %>%
               dplyr::mutate(ancestry = ifelse(ancestry == "parv", "parviglumis", ancestry)),
             aes(y = estimate),
             pch = 18, size = 2) +
  # and errorbars for 95% CI around that mean
  # based on bootstrap
  geom_errorbar(
    data = ci_r5 %>%
      dplyr::mutate(ancestry = ifelse(ancestry == "parv", "parviglumis", ancestry)),
    aes(ymin = low,
        ymax = high),
    width = .5) +
  xlab("Recombination rate quintile (cM/Mbp)") +
  ylab("Proportion ancestry genomewide (HMM)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None") +
  facet_grid(paste(ancestry, "\nancestry") ~ paste("sympatric", zea))
# p_r5_local_anc
ggsave(file = png_r5_local_anc,
       plot = p_r5_local_anc,
       device = "png",
       width = 7.5, height = 7, 
       units = "in", dpi = 300)

# coding density local ancestry plot
# recombination rate local ancestry plots:
p_cd5_local_anc <- anc_by_wind_and_zea %>%
  dplyr::mutate(ancestry = ifelse(ancestry == "parv", "parviglumis", ancestry)) %>%
  ggplot(., aes(x = bin_cd5, y = p, group = paste(zea, ancestry))) +
  # first plot original point estimates for ind. ancestry
  geom_point(aes(shape = zea,
                 color = ancestry),
             alpha = 0.75,
             position = position_jitter(0.2)) +
  scale_color_manual(values = col_maize_mex_parv) +
  # then add mean for that group
  geom_point(data = ci_cd5 %>%
               dplyr::mutate(ancestry = ifelse(ancestry == "parv", "parviglumis", ancestry)),
             aes(y = estimate),
             pch = 18, size = 2) +
  # and errorbars for 95% CI around that mean
  # based on bootstrap
  geom_errorbar(
    data = ci_cd5 %>%
      dplyr::mutate(ancestry = ifelse(ancestry == "parv", "parviglumis", ancestry)),
    aes(ymin = low,
        ymax = high),
    width = .5) +
  xlab("Gene density (coding bp/cM)") +
  ylab("Proportion ancestry genomewide (HMM)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None") +
  facet_grid(paste(ancestry, "\nancestry") ~ paste("sympatric", zea))

# p_cd5_local_anc
ggsave(file = png_cd5_local_anc,
       plot = p_cd5_local_anc,
       device = "png",
       width = 7.5, height = 7, 
       units = "in", dpi = 300)


p_r_cont_local_anc <- anc_by_wind_and_zea %>%
  dplyr::mutate(ancestry = ifelse(ancestry == "parv", "parviglumis", ancestry)) %>%
  dplyr::filter(ancestry != zea) %>%
  dplyr::mutate(label_introgression = factor(paste(ancestry, "in\n", "sympatric", zea),
                                             ordered = T,
                                             levels = paste(c("mexicana", "maize", "parviglumis", "parviglumis"), 
                                                            "in\n", "sympatric", 
                                                            c("maize", "mexicana", "maize", "mexicana")))) %>%
  ggplot(., aes(x = log10(cM_Mb), y = p, group = label_introgression,
                shape = zea, color = ancestry)) +
  geom_point(alpha = 0.75) +
  scale_color_manual(values = col_maize_mex_parv) +
  xlab("Recombination rate log10(cM/Mbp)") +
  ylab("Introgressed ancestry proportion (HMM)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None") +
  facet_wrap(~label_introgression) + 
  geom_smooth(method = "lm", color = "black")
# p_r_cont_local_anc
ggsave(file = png_r_cont_local_anc,
       plot = p_r_cont_local_anc,
       device = "png",
       width = 7.5, height = 7, 
       units = "in", dpi = 300)
ggsave(file = png_r_cont_local_anc_lzw,
       plot = p_r_cont_local_anc,
       device = "tiff",
       width = 7.5, height = 7, 
       units = "in", dpi = 300,
       compression = "lzw", type = "cairo")


# make tidy table of spearman's rank correlation (and bootstrap ci) results
rho <- ci_boot %>%
  filter(feature %in% c("cor.spearman.r", "cor.spearman.cd")) %>%
  mutate(feature = ifelse(feature == "cor.spearman.r", "recombination rate (cM/Mb)",
                          "gene density (coding bp/cM)"),
         method = "ancestry_hmm",
         group = paste("sympatric", zea, sep = " "),
         resolution = "1cM windows") %>%
  dplyr::mutate(ancestry = ifelse(ancestry == "parv", "parviglumis", ancestry)) %>%
  dplyr::select(group, ancestry, feature, estimate, low, high) %>%
  rename(`Spearman's rho` = estimate, `2.5%` = low, `97.5%` = high) %>%
  arrange(., desc(feature), desc(group), ancestry)

# print table to file for estimates of Spearman's rank correlation
print(xtable(rho, 
             digits = c(1, 0, 0, 0, 3, 3, 3),
             type = "latex", 
             latex.environments = NULL),
      include.rownames = F,
      file = file_spearmans_rho)
print(xtable(rho, 
             digits = c(1, 0, 0, 0, 3, 3, 3),
             type = "latex", 
             latex.environments = NULL),
      include.rownames = F,
      file = file_spearmans_rho_tex)


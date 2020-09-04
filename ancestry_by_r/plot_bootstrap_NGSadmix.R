#!/usr/bin/env Rscript

# This script plots the the bootstrap results for
# ancestry across recombination quintiles
library(dplyr)
library(tidyr)
library(ggplot2)

# load variables from Snakefile
windows_file = snakemake@input[["windows"]] # feature information for 1cM windows
# windows_file = "ancestry_by_r/results/map_pos_1cM_windows.txt"
windows <- read.table(windows_file, header = T, stringsAsFactors = F, sep = "\t")
load(snakemake@input[["r5"]]) # NGSAdmix results by recombination quintile
# load("ancestry_by_r/results/bootstrap_1cM/HILO_MAIZE55/r5_K2.Rdata")
load(snakemake@input[["cd5"]]) # NGSAdmix results by coding bp/cM quintile
# load("ancestry_by_r/results/bootstrap_1cM/HILO_MAIZE55/cd5_K2.Rdata")
source(snakemake@input[["colors"]]) # plotting colors
# source("colors.R")
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
png_facet_r5 = snakemake@output[["png_facet_r"]]
# png_facet_r5 = "ancestry_by_r/plots/K2_by_r_bootstrap_lm_elevation_facet_r.png"
png_color_elev_r5 = snakemake@output[["png_color_elev_r"]]
# png_color_elev_r5 = "ancestry_by_r/plots/K2_by_r_bootstrap_lm_elevation_color_elev.png"


# make plots
# mexicana ancestry for K=2 in sympatric maize and mexicana:
p_r5_symp <- r5$anc_group_confidence_intervals %>%
  filter(ancestry == "mexicana_ancestry") %>%
  filter(symp_allo == "sympatric") %>%
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
  # and errorbars for 90% CI around that mean
  # based on bootstrap with NGSAdmix
  geom_errorbar(aes(ymin = low,
                    ymax = high),
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
p_cd5_symp <- cd5$anc_group_confidence_intervals %>%
  filter(ancestry == "mexicana_ancestry") %>%
  filter(symp_allo == "sympatric") %>%
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
  # and errorbars for 90% CI around that mean
  # based on bootstrap with NGSAdmix
  geom_errorbar(aes(ymin = low,
                    ymax = high),
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
p_r5_symp_allo <- r5$anc_group_confidence_intervals %>%
  filter(ancestry == "mexicana_ancestry") %>%
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
  # and errorbars for 90% CI around that mean
  # based on bootstrap with NGSAdmix
  geom_errorbar(aes(ymin = low,
                    ymax = high),
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

p_cd5_symp_allo <- cd5$anc_group_confidence_intervals %>%
  filter(ancestry == "mexicana_ancestry") %>%
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
  # and errorbars for 90% CI around that mean
  # based on bootstrap with NGSAdmix
  geom_errorbar(aes(ymin = low,
                    ymax = high),
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
  filter(., symp_allo == "sympatric") %>%
  filter(., ancestry == "mexicana") %>%
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
#p_facet_r
ggsave(file = png_facet_r5,
       plot = p_facet_r5,
       device = "png",
       width = 7, height = 4, 
       units = "in", dpi = 300)


# different way of visualizing the same pattern:
p_color_elev_r5 <- r5$anc_ind %>%
  filter(., symp_allo == "sympatric") %>%
  filter(., ancestry == "mexicana") %>%
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

#!/usr/bin/env Rscript

# This script plots the the bootstrap results for
# ancestry across recombination quintiles
library(dplyr)
library(tidyr)
library(ggplot2)

# load variables from Snakefile
load(snakemake@params[["meta"]])
# load("samples/HILO_MAIZE55_meta.RData")
source(snakemake@params[["colors"]])
# source("colors.R")
png_symp = snakemake@output[["png_symp"]]
# png_symp = "ancestry_by_r/plots/K2_by_r_bootstrap_sympatric_only.png"
png_symp_allo = snakemake@output[["png_symp_allo"]]
# png_symp_allo = "ancestry_by_r/plots/K2_by_r_bootstrap_sympatric_and_allopatric.png"
png_facet_r = snakemake@output[["png_facet_r"]]
# png_facet_r = "ancestry_by_r/plots/K2_by_r_bootstrap_lm_elevation_facet_r.png"
png_color_elev = snakemake@output[["png_color_elev"]]
# png_color_elev = "ancestry_by_r/plots/K2_by_r_bootstrap_lm_elevation_color_elev.png"

PREFIX = snakemake@params[["prefix_all"]]
# PREFIX = "HILO_MAIZE55"
K = snakemake@params[["k"]]
# K = 2
windows_file = snakemake@params[["windows"]]
# windows_file = "ancestry_by_r/results/map_pos_1cM_windows.txt"

ancestries <- c("maize", "mexicana", "parviglumis")[1:K] # ancestries

alpha = 0.1 # use 90% confidence intervals for bootstrap

# get the recombination rate quintile ranges
r5 <- read.table(windows_file, header = T, stringsAsFactors = F, sep = "\t") %>%
  dplyr::select(., quintile_r5, bin_r5) %>%
  rename(r5_quintile = quintile_r5, r5_bin = bin_r5) %>%
  filter(., !duplicated(r5_bin)) %>%
  arrange(., r5_quintile) %>% # order 1-5
  mutate(r5_bin = factor(r5_bin, ordered = T, levels = r5_bin))


# get bootstrap estimates for proportion ancestry K=2
anc_boot <- do.call(rbind,
                    lapply(0:100, function(BOOT) do.call(rbind,
                                                         lapply(1:5, function(r)
                                                         read.table(paste0("ancestry_by_r/results/bootstrap_1cM/",
                                                         PREFIX, "/r5_", r, "/K", K, "/boot", BOOT, ".anc"),
                                                         header = T, sep = "\t", stringsAsFactors = F) %>%
                                                         mutate(bootstrap = BOOT) %>%
                                                         mutate(r5_quintile = r))))) %>%
  left_join(., meta, by = "ID") %>%
  left_join(., r5, by = "r5_quintile") # label recombination quintile with rates for that bin

# point estimate individual ancestry for each sample
anc_ind <- anc_boot %>%
  filter(bootstrap == 0) %>% # the original sample is called 'bootstrap 0'
  gather(., key = "ancestry", value = "p", ancestries[1:K])

# mean of bootstrap for groups defined by symp/allo, recombination rate and zea subspecies
anc_boot_mean <- anc_boot %>%
  group_by(group, zea, symp_allo, bootstrap, r5_bin) %>%
  summarise(mexicana_ancestry = mean(mexicana),
            maize_ancestry = mean(maize))


# group means (by symp/allo and zea subspecies, bootstrap = 0 is all the original data)
anc_group_estimate <- anc_boot_mean %>%
  filter(bootstrap == 0) %>%
  gather(., key = "ancestry", value = "p",
         paste(ancestries[1:K], "ancestry", sep = "_"))


# get percentiles from bootstrap
anc_boot_perc <- anc_boot_mean %>%
  filter(bootstrap != 0) %>% # this is the original sample
  gather(., key = "ancestry", value = "p",
         paste(ancestries[1:K], "ancestry", sep = "_")) %>%
  group_by(ancestry, group, r5_bin) %>%
  summarise(low_boot = quantile(p, alpha/2), # low and high bounds
            high_boot = quantile(p, 1-alpha/2), # of 90% conf. interval
            median_boot = median(p),
            mean_boot = mean(p))

# calculate basic bootstrap (i.e. 'pivot') confidence intervals together (to plot later as range)
anc_group_confidence_intervals <- anc_group_estimate %>%
  left_join(., anc_boot_perc, by = c("ancestry", "group", "r5_bin")) %>%
  mutate(low = 2*p - high_boot, # p is the estimate from the original sample
         high = 2*p - low_boot)

# mexicana ancestry for K=2 in sympatric maize and mexicana:
p_k2_symp <- anc_group_confidence_intervals %>%
  filter(ancestry == "mexicana_ancestry") %>%
  filter(symp_allo == "sympatric") %>%
  ggplot(aes(x = r5_bin, y = p, group = zea)) +
  # first plot original point estimates for ind. ancestry
  geom_point(data = filter(anc_ind, zea != "parviglumis" & ancestry == "mexicana"),
             aes(x = r5_bin,
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
#p_k2_symp
ggsave(file = png_symp,
       plot = p_k2_symp,
       device = "png",
       width = 5, height = 4, 
       units = "in", dpi = 300)


# mexicana ancestry for K=2 in sympatric and allopatric maize/mex
p_k2_symp_allo <- anc_group_confidence_intervals %>%
  filter(ancestry == "mexicana_ancestry") %>%
  ggplot(aes(x = r5_bin, y = p, group = zea)) +
  # first plot original point estimates for ind. ancestry
  geom_point(data = filter(anc_ind, zea != "parviglumis" & ancestry == "mexicana"),
             aes(x = r5_bin,
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
#p_k2_symp_allo
ggsave(file = png_symp_allo,
       plot = p_k2_symp_allo,
       device = "png",
       width = 7, height = 7, 
       units = "in", dpi = 300)


# plot ancestry across recombination bins and elevation:
p_facet_r <- anc_ind %>%
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
  facet_wrap(~r5_bin) +
  xlab("Elevation (m)") +
  ylab("Proportion mexicana ancestry") +
  guides(color = guide_legend("Subspecies"),
         shape = guide_legend("Subspecies"),
         fill = guide_legend("Subspecies"))
#p_facet_r
ggsave(file = png_facet_r,
       plot = p_facet_r,
       device = "png",
       width = 7, height = 4, 
       units = "in", dpi = 300)


# different way of visualizing the same pattern:
p_color_elev <- anc_ind %>%
  filter(., symp_allo == "sympatric") %>%
  filter(., ancestry == "mexicana") %>%
  ggplot(., aes(x = r5_bin, y = p,
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
#p_color_elev
ggsave(file = png_color_elev,
       plot = p_color_elev,
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

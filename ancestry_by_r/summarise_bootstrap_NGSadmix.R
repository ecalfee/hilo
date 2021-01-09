#!/usr/bin/env Rscript

# This script summarises the bootstrap results for
# ancestry across genomic quintiles
library(dplyr)
library(tidyr)

# load variables from Snakefile
load(snakemake@params[["meta"]])
# load("samples/HILO_MAIZE55_meta.RData")
FEATURE = snakemake@params[["feature"]]
# FEATURE = "r" # vs. "cd"
PREFIX = snakemake@params[["prefix_all"]]
# PREFIX = "HILO_MAIZE55"
K = snakemake@params[["k"]]
# K = 2
windows_file = snakemake@params[["windows"]]
# windows_file = "ancestry_by_r/results/map_pos_1cM_windows.txt"
alpha = snakemake@params[["alpha"]]
# alpha = 0.05 # use 95% confidence intervals for bootstrap
rdata_out = snakemake@output[["rdata"]]
# rdata_out = "ancestry_by_r/results/bootstrap_1cM/HILO_MAIZE55/r5_K2.Rdata"

ancestries <- c("maize", "mexicana", "parviglumis")[1:K] # ancestries

bin_col <- paste0("bin_", FEATURE, "5")
quintile_col <- paste0("quintile_", FEATURE, "5")

# get the genomic quintile ranges
q <- read.table(windows_file, header = T, stringsAsFactors = F, sep = "\t") %>%
  dplyr::select(., quintile_col, bin_col) %>%
  rename(quintile = quintile_col, bin = bin_col) %>%
  filter(., !duplicated(bin)) %>%
  arrange(., quintile) %>% # order 1-5
  mutate(bin = factor(bin, ordered = T, levels = bin))


# get bootstrap estimates for proportion ancestry K=2
anc_boot <- do.call(rbind,
                    lapply(0:100, function(BOOT) do.call(rbind,
                                                         lapply(1:5, function(r)
                                                         read.table(paste0("ancestry_by_r/results/bootstrap_1cM/",
                                                         PREFIX, "/", FEATURE, "5_", r, "/K", K, "/boot", BOOT, ".anc"),
                                                         header = T, sep = "\t", stringsAsFactors = F) %>%
                                                         mutate(bootstrap = BOOT,
                                                         quintile = r,
                                                         feature = FEATURE))))) %>%
  left_join(., meta, by = "ID") %>%
  left_join(., q, by = "quintile") # label genomic quintile with rates for that bin

# point estimate individual ancestry for each sample
anc_ind <- anc_boot %>%
  filter(bootstrap == 0) %>% # the original sample is called 'bootstrap 0'
  gather(., key = "ancestry", value = "p", ancestries[1:K])

# mean of bootstrap for groups defined by symp/allo, recombination rate and zea subspecies
anc_boot_mean <- anc_boot %>%
  group_by(group, zea, symp_allo, bootstrap, bin, quintile, feature) %>%
  summarise(mexicana_ancestry = mean(mexicana),
            maize_ancestry = mean(maize))


# group means (by symp/allo and zea subspecies, bootstrap = 0 is all the original data)
anc_group_estimate <- anc_boot_mean %>%
  filter(bootstrap == 0) %>%
  gather(., key = "ancestry", value = "p",
         paste(ancestries[1:K], "ancestry", sep = "_"))


# calculate percentiles confidence intervals from bootstrap
anc_boot_perc <- anc_boot_mean %>%
  filter(bootstrap != 0) %>% # boot = 0 is the original sample
  gather(., key = "ancestry", value = "p",
         paste(ancestries[1:K], "ancestry", sep = "_")) %>%
  group_by(ancestry, group, bin, quintile, feature) %>%
  summarise(low_boot = quantile(p, alpha/2), # low and high bounds
            high_boot = quantile(p, 1-alpha/2), # of conf. interval
            median_boot = median(p),
            mean_boot = mean(p)) %>%
  left_join(anc_group_estimate, ., by = c("ancestry", "group", "bin", "quintile", "feature"))

# calculate Spearman's rank correlation mexicana ancestry ~ feature (feature = r/cd)
# and 95% percentile confidence intervals from the bootstrap results
# first calculate mean across individuals in each group (sympatric maize, allopatric mexicana..)..
# (so it's just 5 points in the rank correlation; like f4s)
spearman = anc_boot_mean %>%
  group_by(group, bootstrap) %>%
  summarise(rho = cor(x = mexicana_ancestry, 
                           y = quintile, 
                           method = "spearman")) %>%
  ungroup() %>%
  group_by(group) %>%
  summarise(rho_estimate = rho[bootstrap == 0], # original sample
            boot_low = quantile(rho[bootstrap != 0], alpha/2), # percentiles from bootstrap (excl. original sample)
            boot_high = quantile(rho[bootstrap != 0], 1 - alpha/2))
  #ggplot(., aes(fill = group, x = group, y = rho)) +
  #geom_violin()

# calculate significance of change in slope across elevation (elevation*feature interaction)

# save
feature_name = paste0(FEATURE, "5")
assign(x = feature_name, 
       value = list(anc_group_estimate = anc_group_estimate,
                    anc_ind = anc_ind,
                    anc_boot_mean = anc_boot_mean,
                    anc_boot_perc = anc_boot_perc,
                    spearman = spearman,
                    alpha = alpha))
save(list = c(feature_name), file = rdata_out)


#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)

# load variables from Snakefile
windows = snakemake@input[["windows"]]
# windows = "results/map_pos_1cM_windows.txt"
prefix = snakemake@params[["prefix"]]
# prefix = "HILO_MAIZE55"
n = as.integer(snakemake@params[["n"]]) # number of bootstraps
# n = 100
feature = snakemake@params[["feature"]]
# feature = "r"
seed = snakemake@params[["seed"]]
# seed = 100

map_pos_1cM <- read.table(windows, header = T, sep = "\t", stringsAsFactors = F)

quintile = paste0("quintile_", feature, 5) # e.g. quintile_r5

# make n bootstrap replicate lists of windows:
set.seed(seed) # set seed for replication
for (r in 1:5){ # 5 recombination rate quintiles
  r_windows <- map_pos_1cM$window[map_pos_1cM[ , quintile] == r]
  # boot0 is the original data/windows:
  write.table(x = r_windows,
              file = paste0("ancestry_by_r/results/bootstrap_1cM/", prefix, "/", feature, "5_", r, "/boot", 0, ".list"),
              quote = F, col.names = F, row.names = F, sep = "\t")
  for (b in 1:n){ # 100 bootstrap samples
    # (same size = same # windows as original data
    # but windows are sampled w/ replacement)
    write.table(x = sample(x = r_windows,
                         size = length(r_windows),
                         replace = T),
                file = paste0("ancestry_by_r/results/bootstrap_1cM/", prefix, "/", feature, "5_", r, "/boot", b, ".list"),
                quote = F, col.names = F, row.names = F, sep = "\t")
  }
}

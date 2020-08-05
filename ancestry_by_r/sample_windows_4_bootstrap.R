#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)

# load variables from Snakefile
windows = snakemake@input[["windows"]]
# windows = "results/map_pos_1cM_windows.txt"
prefix_all = snakemake@params[["prefix_all"]]
# prefix_all = "HILO_MAIZE55"
n = as.integer(snakemake@params[["n"]]) # number of bootstraps
# n = 100
map_pos_1cM <- read.table(windows, header = T, sep = "\t", stringsAsFactors = F)

# I may want to exclude windows with < 5 SNPs or something similar. First get SNPs in windows.
# r_windows_excl <- c("W382") # e.g. W382 has no SNPs
r_windows_excl <- character() # empty character vector

# make n bootstrap replicate lists of windows:
set.seed(100) # set seed for replication
for (r in 1:5){ # 5 recombination rate quintiles
  r_windows <- map_pos_1cM$window[map_pos_1cM$bin_r5 == levels(map_pos_1cM$bin_r5)[r]]
  r_windows_incl <- r_windows[!(r_windows %in% r_windows_excl)]
  # boot0 is the original data/windows:
  write.table(x = r_windows_incl,
              file = paste0("ancestry_by_r/results/bootstrap_1cM/", prefix_all, "/r5_recomb", r, "/boot", 0, ".list"),
              quote = F, col.names = F, row.names = F, sep = "\t")
  for (b in 1:n){ # 100 bootstrap samples
    # (same size = same # windows as original data
    # but windows are sampled w/ replacement)
    write.table(x = sample(x = r_windows_incl,
                         size = length(r_windows_incl),
                         replace = T),
                file = paste0("ancestry_by_r/results/bootstrap_1cM/", prefix_all, "/r5_recomb", r, "/boot", b, ".list"),
                quote = F, col.names = F, row.names = F, sep = "\t")
  }
}

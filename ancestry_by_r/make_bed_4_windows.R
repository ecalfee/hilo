#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)

# short script writes a separate bed file for each 1cM window

# load variables from Snakefile
windows_file = snakemake@input[["windows"]]
# windows = "results/map_pos_1cM_windows.txt"

out_dir = "ancestry_by_r/results/BED_1cM/"
# out_dir = "results/BED_1cM/"

windows = read.table(windows_file, sep = "\t", header = T, stringsAsFactors = F)

for (w in windows$window){
  filter(windows, window == w) %>%
    write.table(., file = paste0(out_dir, w, ".bed"),
                sep = "\t", col.names = T, row.names = F, quote = F)
}

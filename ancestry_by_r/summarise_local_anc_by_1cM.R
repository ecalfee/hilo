#!/usr/bin/env Rscript
library(dplyr)

# summarises mean local ancestry per window 
#(local ancestry tracts are already mapped onto windows using bedtools)

# load variables from Snakefile
anc_bed_file = snakemake@input[["bed"]]
# anc_bed_file = "ancestry_by_r/results/local_anc_1cM/HILO_MAIZE55/Ne10000_yesBoot/pop365.bed"
windows_file = snakemake@input[["windows"]]
# windows_file = "ancestry_by_r/results/map_pos_1cM_windows.txt"
pop = snakemake@params[["pop"]]
# pop = "pop365"
output_file = snakemake@output[["anc"]]
# output_file = "ancestry_by_r/results/local_anc_1cM/HILO_MAIZE55/Ne10000_yesBoot/pop365.anc.wind"

# windows
windows = read.table(windows_file, sep= "\t", header = T, stringsAsFactors = F) %>%
  dplyr::select(chr, start, end, window)

# summarise local ancestry by window
anc = read.table(anc_bed_file, sep = "\t", header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("chr", "start", "end", "pos", "anc", "window")) %>%
  dplyr::mutate(length = end - start) %>%
  dplyr::group_by(window) %>%
  dplyr::summarise(anc = sum(anc*length)/sum(length)) %>%
  dplyr::left_join(windows, anc, by = "window") %>%
  dplyr::mutate(pop = pop) %>%
  dplyr::select(window, pop, anc)

# write output
write.table(anc, output_file, 
            col.names = F, row.names = F, sep = "\t", quote = F) 


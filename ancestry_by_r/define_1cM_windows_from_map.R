#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)

# load variables from Snakefile
file_out = snakemake@output["windows"]
# file_out = "results/map_pos_1cM_windows.txt"
rmap_file = snakemake@input["rmap"]
# rmap_file = "../data/linkage_map/ogut_fifthcM_map_agpv4_EXTENDED.txt"

# load the extended map (goes to ends of chr's)
rmap_ext = read.table(rmap_file, stringsAsFactors = F, header = T)


# get 1cM windows and their mean recombination rates
map_pos_1cM_v0 <- do.call(rbind,
                       lapply(1:10, # for each chromosome
                              function(x) data.frame(chr = x,
                                                     # get spaced cM positions each 0.2cM
                                                     # starting with the rounded up cM of bp=1
                                                     pos_cM = seq(ceiling(5*rmap_ext$pos_cM[first(which(rmap_ext$chr == x))])/5,
                                                                  rmap_ext$pos_cM[last(which(rmap_ext$chr == x))],
                                                                  by = 1)) %>%
                                # get position in bp using linear approximation
                                mutate(pos_bp = round(approx(x = rmap_ext$pos_cM[rmap_ext$chr == x],
                                                             y = rmap_ext$pos_bp[rmap_ext$chr == x],
                                                             xout = pos_cM,
                                                             method = "linear")$y, 0)) %>%
                                mutate(start = pos_bp - 1) %>% # bed coordinates start at 0
                                mutate(end = c(.$start[-1], NA)) %>% # make non-overlapping by ending 1bp short of next start
                                mutate(width_bp = end - start) %>%
                                mutate(cM_Mb = 1/(width_bp*10^-6)) %>%
                                filter(!is.na(end)) # last position just marks end of the last interval, not it's own start
                       )) %>%
  #name intervals, e.g. w9
  mutate(window = paste0("W", 1:nrow(.)))


# ooops, top 5th of cM windows is a lot smaller than top 5th of genome
# I want to define quantiles from physical space covered, not proportion windows

# get quantile bin for recombination rate
sorted_1cM_bins <- map_pos_1cM_v0 %>%
  arrange(cM_Mb) %>%
  mutate(cum_pos_bp = cumsum(width_bp)) %>%
  mutate(cum_percentile_bp = cum_pos_bp/max(cum_pos_bp))

breaks_5r_1cM <- sapply(seq(0, 1, by = .2), function(x) sorted_1cM_bins$cM_Mb[first(which(sorted_1cM_bins$cum_percentile_bp >= x))])

# label bins in data by their quintile
map_pos_1cM <- map_pos_1cM_v0 %>%
  mutate(bin_r5 = cut(cM_Mb,
                      breaks = breaks_5r_1cM,
                      right = T,
                      include.lowest = T)) %>%
  dplyr::select(., c("chr", "start", "end", "window", "cM_Mb", "bin_r5", "pos_cM"))


# print a list of all 1cM windows and their recombination rate quintile
write.table(map_pos_1cM, file_out,
            quote = F, col.names = T, row.names = F, sep = "\t")

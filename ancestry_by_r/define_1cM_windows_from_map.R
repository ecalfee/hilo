#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(bedr)

# load variables from Snakefile
file_out = snakemake@output[["windows"]]
# file_out = "ancestry_by_r/results/map_pos_1cM_windows.txt"
rmap_file = snakemake@input[["rmap"]]
# rmap_file = "linkage_map/results/ogut_2015_rmap_v2_to_v4_EXTENDED.txt"
cds_file = snakemake@input[["cds"]] # coding regions
# cds_file = "data/refMaize/geneAnnotations/Zea_mays.B73_RefGen_v4.41.chr.CDS.bed"
genome_file = snakemake@input[["genome"]] # genome order
# genome_file = "data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"

# load the extended map (goes to ends of chr's)
rmap_ext = read.table(rmap_file, stringsAsFactors = F, header = T)


# get 1cM windows and their mean recombination rates
map_pos_1cM_v0 <- do.call(rbind,
                       lapply(1:10, # for each chromosome
                              function(x) data.frame(chr = x,
                                                     # get spaced cM positions at 1cM spacing
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
                      include.lowest = T),
         quintile_r5 = as.numeric(factor(bin_r5))) %>% # label quintiles 1-5
  dplyr::select(., c("chr", "start", "end", "window", "cM_Mb", "bin_r5", "quintile_r5", "pos_cM")) %>%
  mutate(chr = as.character(chr))

# print out quintiles
print("Recombination rate quintiles")
print(with(map_pos_1cM, unique(paste(quintile_r5, bin_r5, sep = "_"))))

# map coding bp annotations ("CDS") onto windows and calculate coding_bp (per 1 cM because that's the window size)
sorted_1cM_cds <- bedr(
  engine = "bedtools",
  input = list(a = map_pos_1cM, b = cds_file),
  method = "coverage",
  params = paste0("-sorted -g ", genome_file),
  check.chr = F
) %>%
  data.table::setnames(c(colnames(map_pos_1cM), "n_CDS_overlap", "coding_bp", "width_bp", "frac_bp_coding")) %>%
  mutate(start = as.numeric(start), # bedr changes everything to character string
         end = as.numeric(end),
         cM_Mb = as.numeric(cM_Mb),
         quintile_r5 = as.numeric(quintile_r5),
         pos_cM = as.numeric(pos_cM),
         coding_bp = as.numeric(coding_bp),
         width_bp = as.numeric(width_bp),
         frac_bp_coding = as.numeric(frac_bp_coding)) %>%
  arrange(coding_bp) %>% # arrange from lowest to highest coding base pairs in 1cM window
  # calculate cumulative bp position
  mutate(cum_pos_bp = cumsum(width_bp),
         cum_percentile_bp = cum_pos_bp/max(cum_pos_bp))

# get quantile bin for coding bp/cM
breaks_cd5_1cM <- sapply(seq(0, 1, by = .2),
                         function(x) sorted_1cM_cds$coding_bp[first(which(sorted_1cM_cds$cum_percentile_bp >= x))])

# label bins in data by their quintile
map_pos_1cM_cds <- sorted_1cM_cds %>%
  mutate(bin_cd5 = cut(coding_bp,
                      breaks = breaks_cd5_1cM,
                      dig.lab = 10, # number of digits to display
                      right = T,
                      include.lowest = T),
         quintile_cd5 = as.numeric(factor(bin_cd5))) %>% # label quintiles 1-5
  arrange(chr, start) %>%
  dplyr::select(., c("chr", "start", "end", "window", "cM_Mb", "bin_r5", "quintile_r5", "pos_cM",
                     "coding_bp", "width_bp", "frac_bp_coding", "bin_cd5", "quintile_cd5"))

# print out quintiles
print("coding bp in 1 cM window quintiles")
print(with(map_pos_1cM_cds, unique(paste(quintile_cd5, bin_cd5, sep = "_"))))

# add quintiles for fraction of bp that are coding within the window
sorted_1cM_frac <- map_pos_1cM_cds %>%
  arrange(frac_bp_coding) %>% # arrange from lowest to highest coding base pairs in 1cM window
  # calculate cumulative bp position
  mutate(cum_pos_bp = cumsum(width_bp),
         cum_percentile_bp = cum_pos_bp/max(cum_pos_bp))

# get quantile bin for coding bp/cM
breaks_frac5_1cM <- sapply(seq(0, 1, by = .2),
                         function(x) sorted_1cM_frac$frac_bp_coding[first(which(sorted_1cM_frac$cum_percentile_bp >= x))])

# label bins in data by their quintile
map_pos_1cM_cds_frac <- sorted_1cM_frac %>%
  mutate(bin_frac5 = cut(frac_bp_coding,
                       breaks = breaks_frac5_1cM,
                       right = T,
                       include.lowest = T),
         quintile_frac5 = as.numeric(factor(bin_frac5))) %>% # label quintiles 1-5
  arrange(chr, start) %>%
  dplyr::select(., c("chr", "start", "end", "window", "cM_Mb", "bin_r5", "quintile_r5", "pos_cM",
                     "coding_bp", "bin_cd5", "quintile_cd5",
                     "width_bp", "frac_bp_coding", "bin_frac5", "quintile_frac5")) %>%
  mutate(chr = as.character(chr))

# print out quintiles
print("percent coding bp out of total bp in 1 cM window quintiles")
print(with(map_pos_1cM_cds_frac, unique(paste(quintile_frac5, bin_frac5, sep = "_"))))


# print a list of all 1cM windows and their recombination rate, coding/cM and frac coding bp quintiles
write.table(map_pos_1cM_cds_frac, file_out,
            quote = F, col.names = T, row.names = F, sep = "\t")

#!/usr/bin/env Rscript

# Author: Erin Calfee 2018. Updated 03/2021
# working directory: hilo/
# this script drops 'bad' SNPs from the 0.2 cM Ogut 2015 maize 
# recombination map where a few such SNPs appear to be out
# of order, due to map or assembly error
library(dplyr)
library(ggplot2)

# load input and output file names from snakemake 
# (or alternatively use commented out lines below to load each file outside of a snakemake pipeline)

# Marker postions from Ogut 2015 Supporting file S3 (v2 coordinates)
# were converted to reference genome v4 coordinates using 'Assembly Converter' 
# https://plants.ensembl.org/Zea_mays/Tools/AssemblyConverter
rmap_v4 = snakemake@input[["rmap_v4"]]
# rmap_v4 = "data/linkage_map/ogut_2015_v2_coord_from_supp_file_S3_onto_v4.bed"

# map file with all markers before conversion
rmap_v2 = snakemake@input[["rmap_v2"]]
# rmap_v2 = "data/linkage_map/ogut_2015_v2_coord_from_supp_file_S3.bed"

# output files
# cleaned map
rmap_clean_out = snakemake@output[["rmap_clean"]]
# rmap_clean_out = "linkage_map/results/ogut_2015_rmap_v2_to_v4_INCLUDED.txt"

# all markers
rmap_all_out = snakemake@output[["rmap_all"]]
# rmap_all_out = "linkage_map/results/ogut_2015_rmap_v2_to_v4_ALL.RData"

# get linkage map data
rmap = read.table(rmap_v4, 
                stringsAsFactors = F, header = F) %>%
  data.table::setnames(c("chr", "bed_start", "pos_bp", "marker", "pos_cM")) %>%
  dplyr::select(-bed_start) %>%
  dplyr::filter(substr(chr, 1, 9) != "B73V4_ctg") %>% # remove markers that now mapped to unassembled contigs
  dplyr::mutate(chr = as.integer(chr))

# find positions that have a chromosome that matches neither the chrom in front nor behind them
i_missChr = which(lead(rmap$chr) != rmap$chr & lag(rmap$chr) != rmap$chr)

# visually inspect context for markers that map to the 'wrong' chromosome
# lapply(i_missChr, function(i) rmap[(i-3):(i+3), ])
# remove these SNPs
rmap_0 = rmap[-i_missChr, ]
# confirmed no more
# length(which(lead(rmap_0$chr) != rmap_0$chr & lag(rmap_0$chr) != rmap_0$chr)) == 0
# difference in chromosome number between current and next mapped position
# table(lead(rmap_0$chr) - rmap_0$chr)
i_missChr2 = which(!((lead(rmap_0$chr) - rmap_0$chr) %in% c(0, 1, NA)))
# lapply(i_missChr2, function(i) rmap_0[(i-10):(i+10), ])
# clearly a small segment of chr7 stuck within chr2 
# and a small segment of chr2 stuck within chr10 
# that I'll remove by hand:
rmap_1 = filter(rmap_0, !(marker %in% paste0("M", c(1786:1788, 7190:7196))))
# table(lead(rmap_1$chr) - rmap_1$chr) # looks good! All markers out of order due to chromosome are now removed

# first find reversals
rmap_2 = rmap_1 %>%
  mutate(chr_end = c(diff(chr), 1) == 1,
         length_bp = c(diff(pos_bp), NA),
         length_cM = c(diff(pos_cM), NA),
         length_bp = ifelse(chr_end, NA, length_bp),
         length_cM = ifelse(chr_end, NA, length_cM),
         marker_number = as.integer(substr(marker, 2, 100)),
         reversed = length_bp < 0 & !is.na(length_bp),
         reversed = marker != "M1" & (reversed | lag(reversed)))
# then find smaller out-of-place scaffolds (may look like reversals?)

plot_reversal = function(map, problem_marker, context = 10) {
  map %>%
  mutate(same_chr = chr == map$chr[map$marker_number == problem_marker]) %>%
  filter(same_chr & marker_number > problem_marker - context & marker_number < problem_marker + context) %>%
  ggplot(., aes(x = pos_cM, y = pos_bp, color = reversed)) +
  geom_point(size = .1) +
  theme_light() +
  geom_text(aes(label = marker)) +
  ggtitle(paste("M", problem_marker))
}
#plot_reversal(rmap_2, 37) # good 
#plot_reversal(rmap_2, 5147) # good
#plot_reversal(rmap_2, 770) # not really a reversal

rmap_2_rm_reversals = filter(rmap_2, !reversed)
# which ones are still problems? how many?
non_rev_markers = rmap_2_rm_reversals$marker_number[c(diff(rmap_2_rm_reversals$pos_bp), 1) < 0 & c(diff(rmap_2_rm_reversals$chr), NA) %in% c(0, NA)]
# for this small set of non-reversed problem markers,
# I visualize and select markers to drop by eye
# for (x in non_rev_markers) plot(plot_reversal(map = rmap_2, problem_marker = x, context = 30))
markers_2_remove = paste0("M", c(6465:6466, 4397:4404, 4237:4241, 3966:3967, 1855:1860, 767:769))

original_markers_v2 = read.table(rmap_v2, header = T) %>%
  rename(pos_cM = cM_pos) %>%
  mutate(pos_bp = NA, # don't keep v2 positions
         marker_status = ifelse(marker %in% rmap_1$marker, "mapped", "unmapped or wrong chr")) %>%
  dplyr::select(chr, pos_bp, pos_cM, marker, marker_status)
#table(original_markers_v2$marker_status)/nrow(original_markers_v2)*100
  
# label all markers by their status (keep, unmapped or wrong chr, removed - reversed etc.)  
rmap_all = rmap_1 %>%
  filter(!(marker %in% markers_2_remove)) %>%
  mutate(chr_end = c(diff(chr), 1) == 1,
         length_bp = c(diff(pos_bp), NA),
         length_cM = c(diff(pos_cM), NA),
         length_bp = ifelse(chr_end, NA, length_bp),
         length_cM = ifelse(chr_end, NA, length_cM),
         reversed = length_bp < 0 & !is.na(length_bp),
         reversed = marker != "M1" & (reversed | lag(reversed)),
         marker_status = ifelse(reversed, "remove - reversed", "keep")) %>%
  bind_rows(., filter(rmap_1, marker %in% markers_2_remove) %>%
              mutate(marker_status = "remove - other")) %>%
  bind_rows(., filter(original_markers_v2, marker_status == "unmapped or wrong chr")) %>%
  mutate(marker_number = as.integer(substr(marker, 2, 100))) %>%
  arrange(marker_number)

# save all the markers
save(list = c("rmap_all", "markers_2_remove"), file = rmap_all_out)

# filter map
rmap_keep <- filter(rmap_all, 
       marker_status == "keep")

# calculate recombination rate stats for filtered map
rmap_keep_stats <- rmap_keep %>%
  mutate(
    chr_end = c(diff(chr), 1) == 1,
    length_bp = c(diff(pos_bp), NA),
    length_cM = round(c(diff(pos_cM), NA), 1), # rounds off very small numerical erros (e.g. .199999999). Markers should be every 0.2cM
    length_bp = ifelse(chr_end, NA, length_bp),
    length_cM = ifelse(chr_end, NA, length_cM),
    cM_Mbp = round(length_cM/(length_bp/10^6), 4))
# summary(rmap_keep$cM_Mbp)
# table(rmap_keep$length_cM) # biggest gap is from the long reversal on chr7
# rmap_keep[rmap_keep$length_cM == max(rmap_keep$length_cM, na.rm = T), ]
# View(rmap_keep %>% arrange(., -cM_Mbp))

# make a linkage map output file for all included markers
write.table(rmap_keep,
            file = rmap_clean_out,
            sep = "\t", 
            col.names = T, 
            row.names = F, 
            quote = F)

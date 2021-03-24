#!/usr/bin/env Rscript
# Author: Erin Calfee 2018
# working directory: hilo/
# this script drops 'bad' SNPs from the 0.2 cM Ogut 2015 maize 
# recombination map where a few such SNPs appear to be out
# of order, due to map or assembly error
library(dplyr)
library(ggplot2)

# Marker postions from Ogut 2015 Supporting file S3 (v2 coordinates)
# were converted to reference genome v4 coordinates using 'Assembly Converter' 
# https://plants.ensembl.org/Zea_mays/Tools/AssemblyConverter
rmap_file = "data/linkage_map/ogut_2015_v2_coord_from_supp_file_S3_onto_v4.bed"

rmap_clean_out = "linkage_map/results/ogut_2015_rmap_v2_to_v4_INCLUDED.txt"

rmap = read.table(rmap_file, 
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
rmap_3 = rmap_1 %>%
  filter(!(marker %in% markers_2_remove)) %>%
  mutate(chr_end = c(diff(chr), 1) == 1,
         length_bp = c(diff(pos_bp), NA),
         length_cM = c(diff(pos_cM), NA),
         length_bp = ifelse(chr_end, NA, length_bp),
         length_cM = ifelse(chr_end, NA, length_cM),
         reversed = length_bp < 0 & !is.na(length_bp),
         reversed = marker != "M1" & (reversed | lag(reversed)),
         marker_status = ifelse(reversed, "remove - reversed", "keep")) %>%
  bind_rows(filter(rmap_1, marker %in% markers_2_remove) %>%
              mutate(marker_status = "remove - other")) %>%
  mutate(marker_number = as.integer(substr(marker, 2, 100))) %>%
  arrange(marker_number)

# add recombination rates to filtered map
rmap_keep <- rmap_3 %>%
  dplyr::filter(marker_status == "keep") %>%
  mutate(
    pos_cM = pos_cM, 
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

# make a linkage map output file for included markers
write.table(dplyr::select(rmap_keep, marker, chr, pos_bp, pos_cM, cM_Mbp, length_bp, length_cM, chr_end),
            rmap_clean_out,
            sep = "\t", 
            col.names = T, 
            row.names = F, 
            quote = F)

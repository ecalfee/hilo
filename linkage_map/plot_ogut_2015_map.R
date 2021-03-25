#!/usr/bin/env Rscript
# Author: Erin Calfee. Last updated 03/2021.

# script to plot Ogut 2015 map and dropped SNPs (see clean_ogut_2015_map.R)
# working directory: hilo/
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)

rmap_file = snakemake@input[["rmap_all"]]
# rmap_file = "linkage_map/results/ogut_2015_rmap_v2_to_v4_ALL.RData"
png_out = snakemake@output[["png"]]
# png_out = "linkage_map/plots/ogut_2015_v2_to_v4_rmap.png"
png_out_lzw = snakemake@output[["png_lzw"]]
# png_out_lzw = "../hilo_manuscript/figures_supp/ogut_2015_v2_to_v4_rmap.tif"

# load map data
load(rmap_file)

colors_status = c("blue", "orange", "red", "yellow")
names(colors_status) = c("keep", "remove - reversed", "remove - other", "unmapped or wrong chr")
plot_marker = function(map, problem_marker, label_markers = T, context = 10) {
  focal_chr = map$chr[map$marker_number == problem_marker]
  p = map %>%
    dplyr::mutate(same_chr = sapply(map$chr, function(x) x == focal_chr)) %>%
    filter(!is.na(pos_bp) & same_chr & marker_number > problem_marker - context & marker_number < problem_marker + context) %>%
    ggplot(., aes(y = pos_cM, x = pos_bp/10^6, color = marker_status)) +
    theme_light() +
    xlab("position cM") +
    ylab("position Mbp") +
    ggtitle(paste0("chr", focal_chr, ": M", problem_marker)) +
    scale_color_manual(values = colors_status[1:3]) +
    labs(color = "Marker Status")
  if (label_markers) return(p + geom_point(size = .1) + geom_text(aes(label = marker)))
  else return(p + geom_point())
}
#plot_marker(map = rmap_all, label_markers = F, problem_marker = 37)
#plot_marker(map = rmap_all, problem_marker = 5414, context = 25) # large chunk
#plot_marker(map = rmap_all, problem_marker = 770, context = 10) 
#plot_marker(map = rmap_all, problem_marker = 3452, context = 10) 
#plot_marker(map = rmap_all, problem_marker = 7204, context = 10) 
#plot_marker(map = rmap_all, problem_marker = 5433, context = 50) 
#plot_marker(map = rmap_all, problem_marker = 1980, context = 15) # Jerry (high)
#plot_marker(map = rmap_all, problem_marker = 1228, context = 15) # Emily (high)
#plot_marker(map = rmap_all, problem_marker = 3455, context = 15) # Sal (avg)

# basic stats on the number of markers dropped
p_barplot = ggplot(rmap_all,
                   aes(x = marker_status,
                       fill = marker_status)) + # histogram of # consecutive markers in a reversed section)
  geom_bar() + 
  coord_flip() +
  theme_light() +
  scale_fill_manual(values = colors_status) +
  ylab("marker count") +
  xlab("") +
  theme(legend.position = "None")

print("Number of markers kept/dropped from original Ogut 2015 v2 set")
print(table(rmap_all$marker_status))
print("Percentage of markers kept/dropped from original Ogut 2015 v2 set")
print(table(rmap_all$marker_status)/nrow(rmap_all)*100)

# how many markers in a row are dropped?
rle_all = rle(rmap_all$marker_status == "keep") # counts sets of kept and dropped markers in a row

p_hist = ggplot(data.frame(lengths = rle_all$lengths,
                           values = rle_all$values) %>% 
                  filter(values == F),
                aes(x = lengths)) + # histogram of # consecutive markers in a reversed section)
  geom_histogram(binwidth = 1) + 
  theme_light() +
  xlab("number of markers dropped in a row")



# from all markers that map to the correct chromosome on v4
# look at any where more than 5 in a row are dropped
rmap_v4 <- filter(rmap_all, marker_status != "unmapped or wrong chr")
rle_v4 = rle(rmap_v4$marker_status == "keep") # counts sets of kept and dropped markers in a row

# plot all places where more than 5 markers were dropped in a row
long_drop_markers = rmap_v4$marker_number[cumsum(rle_v4$lengths)[rle_v4$lengths > 5 & !rle_v4$values]] # find sets of consecutives FALSEs with length > 5
list_2_plot_long = lapply(long_drop_markers, function(x)
  ggplotGrob(plot_marker(map = rmap_v4, label_marker = F, problem_marker = x, context = 50) +
               theme(legend.position = "None")))
list_2_plot_long[[length(list_2_plot_long) + 1]] = cowplot::get_legend(plot_marker(map = rmap_all, label_marker = F, problem_marker = long_drop_markers[2]))
p_long <- grid.arrange(grobs = list_2_plot_long,
                       nrow = 4, ncol = 2)
# p_long

# plot all markers removed by hand (not necessarily >5 in a row)?
markers_dropped_by_hand <- as.integer(substr(markers_2_remove, 2, 100))
# find the last marker in each contiguous region
regions_dropped_by_hand <- markers_dropped_by_hand[c(which(diff(markers_dropped_by_hand) != 1), length(markers_dropped_by_hand))]
list_2_plot_other = lapply(regions_dropped_by_hand, function(x)
  ggplotGrob(plot_marker(map = rmap_all, label_marker = F, problem_marker = x, context = 30) +
               theme(legend.position = "None")))
list_2_plot_other[[length(list_2_plot_other) + 1]] = 
  cowplot::get_legend(plot_marker(map = rmap_all, label_marker = F, problem_marker = long_drop_markers[2]) +
                        theme(legend.position = "bottom"))
p_other <- grid.arrange(grobs = list_2_plot_other,
                        layout_matrix = rbind(c(1,2),
                                              c(3,4),
                                              c(5,6),
                                              c(7,7)),
                        heights = c(1,1,1,.2),
                        widths = c(1,1))

# filtering doesn't change the general recombination map/pattern
# zoom in on a pos_cM to pos_bp map. from the general shape of these maps it makes more sense
# to extend the tail recombination rate to markers beyond the map boundaries (on end of chroms)
# than it does to use a chromosome-average for these markers because rates are faster on the ends
# than near the centromere, so this should be a better proxy
p_genome = rmap_all %>%
  filter(!is.na(pos_bp)) %>%
  ggplot(., aes(pos_bp/10^6, pos_cM, color = marker_status)) +
  geom_point(size = 0.1, shape = 22) +
  facet_wrap(.~chr, nrow = 2) +
  theme_light() +
  ylab("cM position") +
  xlab("Mbp position") +
  scale_color_manual(values = colors_status) +
  theme(legend.position = "bottom")

# make combined plot
list_grobs_combined <- list_2_plot_other
# length(list_2_plot_other)
list_grobs_combined[[length(list_2_plot_other) + 1]] <- ggplotGrob(p_genome + theme(legend.position = "None"))
list_grobs_combined[[length(list_2_plot_other) + 2]] <- ggplotGrob(p_barplot)
list_grobs_combined[[length(list_2_plot_other) + 3]] <- ggplotGrob(p_hist)
list_grobs_combined[[length(list_2_plot_other) + 4]] <- textGrob(label = "A", 
                                                                 just = "top",
                                                                 x = unit(0.5, "lines"), 
                                                                 y = unit(0, "lines"))
list_grobs_combined[[length(list_2_plot_other) + 5]] <- textGrob(label = "B", 
                                                                 just = "top",
                                                                 x = unit(0.5, "lines"), 
                                                                 y = unit(0, "lines"))
list_grobs_combined[[length(list_2_plot_other) + 6]] <- textGrob(label = "C", 
                                                                 just = "top",
                                                                 x = unit(0.5, "lines"), 
                                                                 y = unit(0, "lines"))
list_grobs_combined[[length(list_2_plot_other) + 7]] <- textGrob(label = "D", 
                                                                 just = "top",
                                                                 x = unit(0.5, "lines"), 
                                                                 y = unit(0, "lines"))

p_combined = grid.arrange(grobs = list_grobs_combined,
                        layout_matrix = rbind(
                          c(11,NA,12,NA),
                          c(9,9,10,10),
                          c(13,NA,14,NA),
                          c(8,8,1,2),
                          c(8,8,3,4),
                          c(8,8,5,6),
                          c(7,7,7,7)),
                        heights = c(.1,1.5,.1,1,1,1,.2),
                        widths = c(1,1,1,1))

# p_combined

# save combined plot as png
ggsave(filename = png_out,
       plot = p_combined, 
       device = "png", 
       height = 7.5, width = 7.5, 
       dpi = 300, units = "in")
# save combined plot as compressed tiff
ggsave(filename = png_out_lzw,
       plot = p_combined, 
       device = "tiff", 
       width = 7.5, height = 7.5, 
       dpi = 300, units = "in", 
       compression = "lzw", type = "cairo")

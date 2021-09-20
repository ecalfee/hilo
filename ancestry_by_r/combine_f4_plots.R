#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

# this script combines 4 f4 plots

# load variables from Snakefile

# mulipanel plot output
png_out = snakemake@output[["png"]]
# png_out = paste0("ancestry_by_r/plots/f4_multi.png")
png_out_lzw = snakemake@output[["png_lzw"]]
# png_out_lzw = "../hilo_manuscript/figures_supp/f4_multi.tif")


# individual plots input
rds_maize_r5 = snakemake@input[["rds_maize_r5"]]
# rds_maize_r5 = "ancestry_by_r/plots/f4_sympatric_maize_pop22_byr5.plot.rds"
rds_maize_cd5 = snakemake@input[["rds_maize_cd5"]]
# rds_maize_cd5 = "ancestry_by_r/plots/f4_sympatric_maize_pop22_bycd5.plot.rds"
rds_mexicana_r5 = snakemake@input[["rds_mexicana_r5"]]
# rds_mexicana_r5 = "ancestry_by_r/plots/f4_sympatric_mexicana_pop22_byr5.plot.rds"
rds_mexicana_cd5 = snakemake@input[["rds_mexicana_cd5"]]
# rds_mexicana_cd5 = "ancestry_by_r/plots/f4_sympatric_mexicana_pop22_bycd5.plot.rds"

# load plots from rds
plot_maize_r5 = readRDS(rds_maize_r5)
plot_maize_cd5 = readRDS(rds_maize_cd5)
plot_mexicana_r5 = readRDS(rds_mexicana_r5)
plot_mexicana_cd5 = readRDS(rds_mexicana_cd5)

# make multi-panel plot
p_multi <- grid.arrange(grobs = list(textGrob(label = "A", 
                                              x = unit(0.5, "lines"), 
                                              y = unit(0.5, "lines")),
                                     textGrob(label = "B", 
                                              x = unit(0.5, "lines"), 
                                              y = unit(0.5, "lines")),
                                     ggplotGrob(plot_maize_r5),
                                     ggplotGrob(plot_maize_cd5),
                                     textGrob(label = "C", 
                                              x = unit(0.5, "lines"), 
                                              y = unit(0.5, "lines")),
                                     textGrob(label = "D", 
                                              x = unit(0.5, "lines"), 
                                              y = unit(0.5, "lines")),
                                     ggplotGrob(plot_mexicana_r5),
                                     ggplotGrob(plot_mexicana_cd5)),
                        layout_matrix = rbind(
                          c(1, 2),
                          c(3, 4),
                          c(5, 6),
                          c(7, 8)),
                        heights = c(0.1, 1, 0.1, 1),
                        widths = c(1, 1))

# p_multi
ggsave(png_out, 
       plot = p_multi, 
       device = "png", 
       width = 7.5, 
       height = 7.5, 
       units = "in",
       dpi = 300)
ggsave(png_out_lzw, 
       plot = p_multi, 
       device = "tiff", 
       width = 7.5, 
       height = 7.5, 
       units = "in",
       compression = "lzw", type = "cairo",
       dpi = 300)

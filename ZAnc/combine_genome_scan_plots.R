#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

# this script combines 4 outlier plots:
# mean ancestry in maize and mexicana sympatric pops
# and slope with elevation

# load variables from Snakefile

# mulipanel plot output
png_out = snakemake@output[["png"]]
# png_out = paste0("ZAnc/plots/Ne10000_yesBoot/multi_maize_mexicana_genome_scan.png")

# individual plots input
rds_maize_mean = snakemake@input[["rds_maize_mean"]]
# rds_maize_mean = paste0("ZAnc/plots/Ne10000_yesBoot/maize.meanAnc.plot.rds")
rds_mexicana_mean = snakemake@input[["rds_mexicana_mean"]]
# rds_mexicana_mean = paste0("ZAnc/plots/Ne10000_yesBoot/mexicana.meanAnc.plot.rds")
rds_maize_elev = snakemake@input[["rds_maize_elev"]]
# rds_maize_elev = paste0("ZAnc/plots/Ne10000_yesBoot/maize.lmElev.plot.rds")
rds_mexicana_elev = snakemake@input[["rds_mexicana_elev"]]
# rds_mexicana_elev = paste0("ZAnc/plots/Ne10000_yesBoot/mexicana.lmElev.plot.rds")

# load plots from rds
plot_maize_mean = readRDS(rds_maize_mean)
plot_mexicana_mean = readRDS(rds_mexicana_mean)
plot_maize_elev = readRDS(rds_maize_elev)
plot_mexicana_elev = readRDS(rds_mexicana_elev)

# make simple plot just to extract legend only for genomewide mean and 5% FDR lines on each plot
plot_for_legend_only <- ggplot(data = data.frame(lines = c("genomewide_mean", "fdr5"), 
                         value = 1)) +
  geom_hline(aes(yintercept = value, color = lines, linetype = lines)) +
  theme_classic() +
  labs(color = "a", linetype = "a") +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  scale_colour_manual(values = c(genomewide_mean = "black",
                                 fdr5 = "#00BFC4"),
                      labels = c("genomewide mean", "5% FDR"),
                      limits = c("genomewide_mean", "fdr5")) +
  scale_linetype_manual(values = c(genomewide_mean = "dashed",
                                 fdr5 = "solid"),
                      labels = c("genomewide mean", "5% FDR"),
                      limits = c("genomewide_mean", "fdr5"))

# make multi-panel plot
p_multi <- grid.arrange(grobs = list(textGrob(label = "A", 
                                              x = unit(0.5, "lines"), 
                                              y = unit(0.5, "lines")),
                                     ggplotGrob(plot_maize_mean + 
                                                  labs(subtitle = "  sympatric maize") +
                                                  guides(color = F) + 
                                                  theme(axis.title.x = element_blank(),
                                                        axis.title.y = element_blank(),
                                                        plot.title = element_blank())),
                                     textGrob(label = "B", 
                                              x = unit(0.5, "lines"), 
                                              y = unit(0.5, "lines")),
                                     ggplotGrob(plot_mexicana_mean + 
                                                  labs(subtitle = "  sympatric mexicana") +
                                                  guides(color = F) + 
                                                  theme(axis.title.x = element_blank(),
                                                        axis.title.y = element_blank(),
                                                        plot.title = element_blank())),
                                     ggplotGrob(plot_maize_elev + 
                                                  labs(subtitle = "  sympatric maize") +
                                                  guides(color = F) + 
                                                  theme(axis.title.x = element_blank(),
                                                        axis.title.y = element_blank(),
                                                        plot.title = element_blank())),
                                     ggplotGrob(plot_mexicana_elev + 
                                                  labs(subtitle = "  sympatric mexicana") +
                                                  guides(color = F) + 
                                                  theme(axis.title.x = element_blank(),
                                                        axis.title.y = element_blank(),
                                                        plot.title = element_blank())),
                                     textGrob(label = "bp position on chromosomes (total length = 2.3Gb)", 
                                              #y = unit(0.5, "lines"),
                                              just = "center",
                                              gp = gpar(fontsize = 11)),
                                     cowplot::get_legend(plot_for_legend_only),
                                     textGrob(label = "mean mexicana ancestry",
                                              rot = 90),
                                     textGrob(label = "slope mexicana ancestry ~ elevation",
                                              rot = 90)
                                     ),
                        layout_matrix = rbind(
                          c(1, 9, NA),
                          c(NA, 9, 2),
                          c(NA, 9, 4),
                          c(3, 10, NA),
                          c(NA, 10, 5),
                          c(NA, 10, 6),
                          c(NA, NA, 8),
                          c(NA, NA, 7)),
                        heights = c(0.15,1,1,0.05,1,1,.1,.25),
                        widths = c(.1, .1, 5))

# p_multi
ggsave(png_out, 
       plot = p_multi, 
       device = "png", 
       width = 7.5, 
       height = 7.5, 
       units = "in",
       dpi = 300)


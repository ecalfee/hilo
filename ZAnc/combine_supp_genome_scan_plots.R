#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

# this script combines 4 outlier plots:
# mean parv/mexicana ancestry in maize and mexicana sympatric pops
# for supplement

# load variables from Snakefile

# mulipanel plot output
png_out = snakemake@output[["png"]]
# png_out = "ZAnc/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_supp_parv_genome_scan.png"
png_out_lzw = snakemake@output[["png_lzw"]]
# png_out_lzw = "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_supp_parv_genome_scan.tif"


# individual plots input
rds_maize_parv = snakemake@input[["rds_maize_parv"]]
# rds_maize_parv = "ZAnc/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_maize_mean_parv_anc.plot.rds"
rds_mexicana_parv = snakemake@input[["rds_mexicana_parv"]]
# rds_mexicana_parv = "ZAnc/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_mexicana_mean_parv_anc.plot.rds"

# other genome coordinates input
genome_file = snakemake@input[["genome"]]
# genome_file = "data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"

# chromosome lengths
genome <- read.table(genome_file, header = F, stringsAsFactors = F,
                     sep = "\t") %>%
  data.table::setnames(c("chr", "length")) %>%
  dplyr::mutate(chr_end = cumsum(length),
                chr_start = c(0, chr_end[1:(nrow(.)-1)]))




# load plots from rds
plot_maize_parv = readRDS(rds_maize_parv)
plot_mexicana_parv = readRDS(rds_mexicana_parv)

# make simple plot just to extract legend only for genomewide mean and FDR lines on each plot
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
p_multi <- grid.arrange(grobs = list(
                                     ggplotGrob(plot_maize_parv + 
                                                  labs(subtitle = "  sympatric maize") +
                                                  guides(color = F) + 
                                                  theme(axis.title.x = element_blank(),
                                                        axis.title.y = element_blank(),
                                                        plot.title = element_blank())),
                                     ggplotGrob(plot_mexicana_parv + 
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
                                     textGrob(label = "mean parviglumis ancestry",
                                              rot = 90)
                                     ),
                        layout_matrix = rbind(
                          c(5, 1),
                          c(5, 2),
                          c(NA, 3),
                          c(NA, 4)),
                        heights = c(2,2,.1,.25),
                        widths = c(.1, 5))

# p_multi
ggsave(png_out, 
       plot = p_multi, 
       device = "png", 
       width = 7.5, 
       height = 4, 
       units = "in",
       dpi = 300)
ggsave(png_out_lzw, 
       plot = p_multi, 
       device = "tiff", 
       width = 7.5, 
       height = 4, 
       units = "in",
       compression = "lzw", type = "cairo",
       dpi = 300)

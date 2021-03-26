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
png_out_lzw = snakemake@output[["png_lzw"]]
# png_out_lzw = "../hilo_manuscript/figures_main/Ne10000_yesBoot_multi_maize_mexicana_genome_scan.tif")


# individual plots input
rds_maize_mean = snakemake@input[["rds_maize_mean"]]
# rds_maize_mean = paste0("ZAnc/plots/Ne10000_yesBoot/maize.meanAnc.plot.rds")
rds_mexicana_mean = snakemake@input[["rds_mexicana_mean"]]
# rds_mexicana_mean = paste0("ZAnc/plots/Ne10000_yesBoot/mexicana.meanAnc.plot.rds")
rds_maize_elev = snakemake@input[["rds_maize_elev"]]
# rds_maize_elev = paste0("ZAnc/plots/Ne10000_yesBoot/maize.lmElev.plot.rds")
rds_mexicana_elev = snakemake@input[["rds_mexicana_elev"]]
# rds_mexicana_elev = paste0("ZAnc/plots/Ne10000_yesBoot/mexicana.lmElev.plot.rds")

# other genome coordinates input
genome_file = snakemake@input[["genome"]]
# genome_file = "data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"
mhl1_QTL_bed = snakemake@input[["mhl1_QTL"]]
# mhl1_QTL_bed = "data/known_QTL/chr9_bin4_mhl1_locus_v4.bed"
inv_file = snakemake@input[["inv"]]
# inv_file = "data/refMaize/inversions/knownInv_v4_coord.txt"
genes_list = snakemake@input[["genes"]]
# "ZAnc/results/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/genes_mapped_to_outliers.txt"
# genes_list = "ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/genes_mapped_to_outliers.txt"


# chromosome lengths
genome <- read.table(genome_file, header = F, stringsAsFactors = F,
                     sep = "\t") %>%
  data.table::setnames(c("chr", "length")) %>%
  dplyr::mutate(chr_end = cumsum(length),
                chr_start = c(0, chr_end[1:(nrow(.)-1)]))

# positions on chr for key loci
inv4m <- read.table(inv_file, header = F, sep = "\t", stringsAsFactors = F) %>%
  data.table::setnames(c("name", "chr", "start", "end", "length")) %>%
  dplyr::filter(name == "inv4m") %>%
  dplyr::mutate(name = "Inv4m") %>%
  dplyr::select(name, chr, start, end)

mhl1 <- read.table(mhl1_QTL_bed, header = F, sep = "\t", stringsAsFactors = F) %>%
  data.table::setnames(c("chr", "start", "end")) %>%
  dplyr::mutate(name = "mhl1")

hpc1 <- read.csv(genes_list, header = T, sep = "\t", stringsAsFactors = F) %>%
  dplyr::rename(name = name_short) %>%
  dplyr::filter(name == "HPC1") %>%
  dplyr::select(name, chr, start, end)

# get (cumulative bp) positions of key loci on the genome
loci = bind_rows(inv4m, mhl1, hpc1) %>%
  left_join(., genome, by = "chr") %>%
  dplyr::mutate(start_cum = start + chr_start,
                end_cum = end + chr_start,
                pos_cum = (start_cum + end_cum)/2,
                length = end - start)


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

# add labels for key loci to maize elevation plot
plot_maize_elev_labelled <- plot_maize_elev + 
             labs(subtitle = "  sympatric maize") +
             guides(color = F) + 
             theme(axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   plot.title = element_blank()) +
  scale_y_continuous(breaks = c(-.5, 0, .5, 1),
                     limits = c(-0.7, 1.45),
                     expand = expansion(mult = c(0, 0),
                                        add = c(0, 0))) +
  geom_text(data = loci, aes(label = name, x = pos_cum), color = "black", y = 1.37, size = 3) +
  geom_segment(data = loci, aes(x = start_cum, xend = end_cum), y = 1.21, yend = 1.21, color = "black") +
  geom_segment(data = loci, aes(x = start_cum, xend = start_cum), y = 1.17, yend = 1.25, color = "black") +
  geom_segment(data = loci, aes(x = end_cum, xend = end_cum), y = 1.17, yend = 1.25, color = "black")
  
# plot_maize_elev_labelled

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
                                     ggplotGrob(plot_maize_elev_labelled),
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
                        heights = c(0.15,1,1,0.05,1.1,1,.1,.25),
                        widths = c(.1, .1, 5))

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

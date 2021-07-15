#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

# this script plots slope from lm ancestry ~ elev
# across the broader mhl1 QTL region and putative inversion

# load variables from Snakefile
zea = "maize"
colors_file = snakemake@input[["colors"]]
# colors_file = "colors.R"
fdr_maize_file = snakemake@input[["fdr_maize"]]
# fdr_maize_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".lmElev.fdr.RData")
fit_maize_file = snakemake@input[["fit_maize"]]
# fit_maize_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".lmElev.fit.RData")
meta_maize_file = snakemake@input[["meta_maize"]]
# meta_maize_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pop.meta.RData")
sites_file = snakemake@input[["sites"]]
# sites_file = "local_ancestry/results/thinnedSNPs/HILO_MAIZE55/whole_genome.var.sites"
png_out = snakemake@output[["png"]]
# png_out = "mhl1_inv/plots/HILO_MAIZE55/Ne10000_yesBoot/mhl1_inv_ancestry.png"
png_out_lzw = snakemake@output[["png_lzw"]]
# png_out_lzw = "../hilo_manuscript/figures_supp/Ne10000_yesBoot_mhl1_inv_ancestry.tif"
mhl1_QTL_bed = snakemake@input[["mhl1_QTL_bed"]]
# mhl1_QTL_bed = "data/known_QTL/chr9_bin4_mhl1_locus_v4.bed"
mhl1_inv_bed = snakemake@input[["mhl1_inv_bed"]]
# mhl1_inv_bed = "mhl1_inv/results/HILO_MAIZE55/Ne10000_yesBoot/mhl1_inv.bed"

# load data
source(colors_file)
load(fdr_maize_file)
load(fit_maize_file)
load(meta_maize_file)


# load site/position information for SNPs with ancestry calls
sites <- read.table(sites_file, header = F, stringsAsFactors = F,
                    sep = "\t") %>%
  data.table::setnames(c("chr", "pos", "major", "minor"))

qtl9 <- read.table(mhl1_QTL_bed,
                   header = F, sep = "\t")
colnames(qtl9) <- c("chr", "start", "end")

p_mhl1 <- bind_cols(sites, fits) %>%
  filter(chr == 9) %>%
  mutate(even_chr = ifelse(chr %% 2 == 0, "even", "odd"),
         zea = zea) %>%
  ggplot(., aes(pos/10^6, envWeights, 
                color = even_chr)) +
  #geom_vline(xintercept = c(inv_mhl1[["start"]], inv_mhl1[["end"]])/10^6, linetype = "dashed", color = "orange") +
  geom_hline(yintercept = filter(FDRs, FDR == 0.05)$thresholds, linetype = "solid", color = "#00BFC4") +
  geom_point(size = .1) +
  geom_hline(yintercept = mean(fits$envWeights), color = "black", linetype = "dashed") +
  xlab("position on chr9 (Mbp)") +
  ylab("slope mexicana ancestry ~ elevation") +
  scale_colour_manual(values = c(odd = "darkgrey", 
                                 even = unname(col_maize_mex_parv[zea]))) + 
  scale_x_continuous(expand = expansion(mult = c(0, 0),
                                        add = c(0, 0)),
                     limits = c(80, 140)) +
  scale_y_continuous(breaks = c(-.5, 0, .5, 1),
                     limits = c(-0.7, 1.45),
                     expand = expansion(mult = c(0, 0),
                                        add = c(0, 0))) +
  theme(legend.position = "none") +
  theme_classic() +
  guides(color = F) +
  geom_text(label = "mhl1 QTL", x = sum(c(qtl9[["start"]], qtl9[["end"]]))/(2*10^6), size = 4,
            color = "black", y = 1.37) +
  geom_segment(x = qtl9[["start"]]/10^6, xend = qtl9[["end"]]/10^6, y = 1.21, yend = 1.21, color = "black") +
  geom_segment(x = qtl9[["start"]]/10^6, xend = qtl9[["start"]]/10^6, y = 1.17, yend = 1.25, color = "black") +
  geom_segment(x = qtl9[["end"]]/10^6, xend = qtl9[["end"]]/10^6, y = 1.17, yend = 1.25, color = "black")
# p_mhl1

# make simple plot just to extract legend only for genomewide mean and 5% FDR lines
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

p <- grid.arrange(grobs = list(ggplotGrob(p_mhl1),
                               cowplot::get_legend(plot_for_legend_only)),
                  nrow = 2,
                  heights = c(1,.1))

# p
ggsave(png_out, 
       plot = p, 
       device = "png", 
       width = 7.5, 
       height = 4, 
       units = "in",
       dpi = 300)

ggsave(png_out_lzw, 
       plot = p, 
       device = "tiff", 
       width = 7.5, 
       height = 4, 
       units = "in",
       dpi = 300,
       compression = "lzw", type = "cairo")


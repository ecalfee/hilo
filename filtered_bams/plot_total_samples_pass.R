#!/usr/bin/env Rscript
# working directory is hilo/
# this script summarises filtered bam metrics and creates a
# plot of samples that pass minimum coverage threshold

library(dplyr)
library(tidyr)
library(ggplot2)
# get colors for plot
# source("colors.R")
source(snakemake@input[["colors"]])
# load sample metadata, including estimated sequence coverage
# load("samples/HILO_MAIZE55_meta.RData")
load(snakemake@input[["meta"]])
# plots out
png = snakemake@output[["png"]]
# png = "filtered_bams/plots/p_seq_counts.png"
png_lzw = snakemake@output[["png_lzw"]]
# png_lzw = "../hilo_manuscript/figures_supp/p_seq_counts.tif"

# plot total included samples per population
p_seq_counts <- meta %>%
  filter(dataset == "HILO") %>%
  mutate(coverage = cut(est_coverage, 
                        breaks = c(0, 0.05, 0.5, 10)
  )) %>%
  ggplot(., aes(alpha = coverage,
                fill = zea,
                x = LOCALITY)) +
  geom_bar() +
  labs(fill = "Zea", x = "Population", y = "# Individuals") +
  scale_alpha_discrete(range = c(.6, 1),
                       name = "Coverage",
                       labels = c(expression("x "<" 0.5"), expression("x ">=" 0.5"))) +
  scale_fill_manual(values = col_maize_mex_parv) +
  facet_grid(zea~.) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
#p_seq_counts
# save plots
ggsave(filename = png,
       plot = p_seq_counts,
       height = 4, width = 5.4, 
       units = "in", device = "png")
ggsave(filename = png_lzw,
       plot = p_seq_counts,
       height = 4, width = 5.4, 
       units = "in", device = "tiff",
       compression = "lzw", type = "cairo")
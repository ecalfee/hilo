#!/usr/bin/env Rscript
# working directory is hilo/
library(GenomicRanges)
library(data.table)
library(dplyr)
library(ggplot2)
# script creates bed file of top X% raisd hits (as contiguous regions) from raw raisd output

# input files
bed_in <- snakemake@input[["hits"]]
# bed_in <- "domestication_scan/results/raisdHits_keep.bed"
overlap_file <- snakemake@input[["overlap"]]
# overlap_file <- "domestication_scan/results/HILO_MAIZE55/Ne10000_yesBoot/raisdOverlap.maize_neg_meanAnc_outliers.perc02.counts"
genome_file <- snakemake@input[["genome"]]
# genome_file <- "data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"
# output files
png_out <- snakemake@output[["png"]]
# png_out <- "domestication_scan/results/HILO_MAIZE55/Ne10000_yesBoot/raisdOverlap.maize_pos_meanAnc_outliers.perc02.png"
png_out_lzw <- snakemake@output[["png_lzw"]]
# png_out_lzw <- "../hilo_manuscript/figures_supp/Ne10000_yesBoot/raisdOverlap.maize_pos_meanAnc_outliers.perc02.tif"


summary_out <- snakemake@output[["summary"]]
# summary_out <- "domestication_scan/results/HILO_MAIZE55/Ne10000_yesBoot/raisdOverlap.maize_pos_meanAnc_outliers.perc02.summary" 

hits <- read.table(bed_in, stringsAsFactors = F, sep = "\t") %>%
  data.table::setnames(c("chr", "start", "end")) %>%
  mutate(length = end - start)

# load in results. First row is the original data, then permutations
overlap <- read.table(overlap_file) %>%
  data.table::setnames(c("n_overlap", "bp_overlap")) %>%
  dplyr::mutate(d = c("original", rep("shuffled", nrow(.) - 1)),
         percent_bp_overlap = 100*bp_overlap/sum(hits$length))

# plot permutation test results
percent_bp_overlap_original = filter(overlap, d == "original")$percent_bp_overlap
p <- overlap %>%
  dplyr::filter(d == "shuffled") %>%
  ggplot(aes(x = percent_bp_overlap)) +
  geom_histogram() +
  geom_vline(xintercept = percent_bp_overlap_original, 
             color = "blue") +
  theme_light() +
  xlab("% of 'domestication' bp overlapping 'introgression deserts'")

ggsave(png_out, 
       plot = p, 
       device = "png", 
       width = 5.4, height = 3.5, units = "in",
       dpi = 300)
ggsave(png_out_lzw, 
       plot = p, 
       device = "tiff", 
       width = 5.4, height = 3.5, units = "in",
       dpi = 300,
       compression = "lzw", type = "cairo")

# summarise results in a table
genome_length_kb <- read.table(genome_file, header = F)$V2 %>%
  sum()/10^3

data.frame(
  percent_bp_overlap = percent_bp_overlap_original,
  p_value = filter(overlap, d == "shuffled") %>%
    summarise(p_value = sum(percent_bp_overlap >= percent_bp_overlap_original)/n()),
  total_hits_kb = sum(hits$length)/1000) %>%
  dplyr::mutate(perc_of_genome_in_hits = 100*total_hits_kb/genome_length_kb) %>%
  write.table(., 
            file = summary_out,
            quote = FALSE, sep = '\t',
            row.names = FALSE, col.names = TRUE)

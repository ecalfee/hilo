#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(xtable)

# this script compares ancestry outliers to known key genes

# load variables from Snakefile
maize_bed = snakemake@input[["maize_bed"]]
# maize_bed = "domestication_scan/results/HILO_MAIZE55/Ne10000_yesBoot/domestication_genes_from_lit.plus20kb.maize.min_mexicana_ancestry.bed"
maize_overlap = snakemake@input[["maize_overlap"]]
# maize_overlap = "domestication_scan/results/HILO_MAIZE55/Ne10000_yesBoot/domestication_genes_from_lit.plus20kb.maize_neg_meanAnc_outliers.perc05.bed"
mexicana_bed = snakemake@input[["mexicana_bed"]]
# mexicana_bed = "domestication_scan/results/HILO_MAIZE55/Ne10000_yesBoot/domestication_genes_from_lit.plus20kb.mexicana.max_mexicana_ancestry.bed"
mexicana_overlap = snakemake@input[["mexicana_overlap"]]
# mexicana_overlap = "domestication_scan/results/HILO_MAIZE55/Ne10000_yesBoot/domestication_genes_from_lit.plus20kb.mexicana_pos_meanAnc_outliers.perc05.bed"

tbl_out = snakemake@output[["tbl"]]
# tbl_out = "domestication_scan/tables/HILO_MAIZE55/Ne10000_yesBoot/domestication_genes.tex"

maize_hits <- read.table(maize_overlap, stringsAsFactors = F, sep = "\t", header = F)$V4
maize <- read.table(maize_bed, sep = "\t", header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("chr", "start", "end", "gene", "min_mex")) %>%
  mutate(min_intro = round(min_mex, 3),
         min_intro_maize = ifelse(gene %in% maize_hits, paste0(min_intro, "*"), min_intro))

mexicana_hits <- read.table(mexicana_overlap, stringsAsFactors = F, sep = "\t", header = F)$V4
mexicana <- read.table(maize_bed, sep = "\t", header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("chr", "start", "end", "gene", "max_mex")) %>%
  mutate(min_intro = round(1 - max_mex, 3),
         min_intro_mex = ifelse(gene %in% mexicana_hits, paste0(min_intro, "*"), min_intro))

combined <- full_join(maize, mexicana, by = c("chr", "start", "end", "gene")) %>%
  arrange(gene) %>%
  mutate(v4_coord = paste0(chr, ":", start + 1, "-", end)) %>%
  dplyr::select(gene, v4_coord, min_intro_maize, min_intro_mex) %>%
  rename(`v4 coordinates` = v4_coord,
         `min introgression\n in maize` = min_intro_maize,
         `min introgression\n in mexicana` = min_intro_mex)

print(xtable(combined,
             type = "latex",
             latex.environments = NULL),
      include.rownames = F,
      file = tbl_out)

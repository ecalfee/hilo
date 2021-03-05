#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)

# this script compares ancestry outliers to known key genes

# load variables from Snakefile
genes_list = snakemake@input[["genes_list"]]
# genes_list = "data/key_genes.csv"
bed_out = snakemake@output[["bed"]]
# bed_out = "domestication_scan/results/domestication_genes_from_lit.bed"

# load list of key genes
genes <- read.csv(genes_list, sep = ",", header = T) %>%
  dplyr::filter(., include) %>%
  tidyr::separate(., v4_coord, c("chr", "start", "end"), remove = F) %>%
  dplyr::mutate(start = as.numeric(start) - 1, # put into bed style 0-based coordinates
                end = as.numeric(end),
                type = "gene")


# isolate just gene coordinates on v4
genes_domestication <- dplyr::filter(genes, category == "domestication") %>%
  dplyr::select(chr, start, end, name_short) %>%
  dplyr::arrange(as.integer(chr), start)


# print bed file of domestication genes only
write.table(genes_domestication, file = bed_out,
            col.names = T, row.names = F, quote = F, sep = "\t")

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
tbl_out_tex = snakemake@output[["tbl_tex"]]
# tbl_out_tex = "../hilo_manuscript/tables/Ne10000_yesBoot_domestication_genes.tex"


maize_hits <- read.table(maize_overlap, stringsAsFactors = F, sep = "\t", header = F)$V4
maize <- read.table(maize_bed, sep = "\t", header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("chr", "start", "end", "gene", "min_mex")) %>%
  mutate(min_intro = round(min_mex, 3),
         min_intro_maize = ifelse(gene %in% maize_hits, paste0(min_intro, "*"), min_intro))

mexicana_hits <- read.table(mexicana_overlap, stringsAsFactors = F, sep = "\t", header = F)$V4
mexicana <- read.table(mexicana_bed, sep = "\t", header = F, stringsAsFactors = F) %>%
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

meta <- rbind(c("zagl1", "ear size" , "cite{Wills:2018_zagl1}"),
              c("gt1", "prolificacy" , "cite{Wills:2013_gt1}"),
              c("ZmSh1-1", "seed shattering" , "cite{Lin:2012_shattering}"),
              c("tb1", "branching" , "cite{Doebley_Stec_Gustus:1995_tb1, Doebley_Stec_Hubbard:1997_tb1, Dong:2019_reg_domestication}"),
              c("zfl2", "cob rank ", "cite{Doebley_Stec:1991, Doebley_Stec:1993, Bomblies_Doebley:2006}"),
              c("pbf1", "storage protein synthesis" , "cite{Wang:1998_pbf}"),
              c("ra2" , "inflorescence architecture" , "cite{Vollbrecht:2005_ramosa}"),
              c("ba1" , "plant architecture" , "cite{Gallavotti:2004_ba1}"),
              c("su1" , "starch biosynthesis" , "cite{Whitt:2002_starch}"),
              c("tga1" , "'naked' grains" , "cite{Dorweiler:1993, Wang:2005_tga1}"),
              c("bt2" , "starch biosynthesis" , "cite{Whitt:2002_starch}"),
              c("ZmSh1-5.1+ZmSh1-5.2" , "seed shattering" , "cite{Lin:2012_shattering}"),
              c("sweet4c" , "sugar transport and seed size" , "cite{Sosso:2015}"),
              c("ae1" , "starch biosynthesis" , "cite{Whitt:2002_starch}"),
              c("ra1" , "inflorescence architecture" , "cite{Vollbrecht:2005_ramosa, Sigmon_Vollbrecht:2010}")
              ) %>%
  as.data.frame(., stringsAsFactors = F) %>%
  data.table::setnames(c("gene", "phenotype", "refs"))

print(xtable(x = left_join(meta, combined, by = "gene"),
             type = "latex",
             latex.environments = NULL),
      include.rownames = F,
      digits = NULL,
      file = tbl_out)

print(xtable(x = left_join(meta, combined, by = "gene"),
             type = "latex",
             latex.environments = NULL),
      include.rownames = F,
      digits = NULL,
      file = tbl_out_tex)

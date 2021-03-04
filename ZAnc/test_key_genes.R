#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(bedr)
library(xtable)

# this script compares ancestry outliers to known key genes

# load variables from Snakefile
genome_file = snakemake@input[["genome"]]
# genome_file = "data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"
Ne = snakemake@params[["Ne"]]
# Ne = 10000
yesno = snakemake@params[["yesno"]]
# yesno = "yes"
prefix = snakemake@params[["prefix"]]
# prefix = "HILO_MAIZE55"
genes_list = snakemake@input[["genes_list"]]
# genes_list = "data/key_genes.csv"
bed_sites_file = snakemake@input[["bed_sites"]]
# bed_sites_file = "local_ancestry/results/thinnedSNPs/HILO_MAIZE55/whole_genome.bed"
txt_out = snakemake@output[["txt"]]
# txt_out = paste0("ZAnc/results/", prefix, "/Ne", Ne, "_", yesno, "Boot/genes_mapped_to_outliers.txt")
tbl_out = snakemake@output[["tbl"]]
# tbl_out = paste0("ZAnc/tables/", prefix, "/Ne", Ne, "_", yesno, "Boot/genes_mapped_to_outliers.tex")

bed_sites <- read.table(bed_sites_file, header = F, 
                        sep = "\t", stringsAsFactors = F) %>%
  data.table::setnames(c("chr", "start", "end", "length")) %>%
  dplyr::mutate(chr = as.character(chr))

# load list of key genes
genes <- read.csv(genes_list, sep = ",", header = T) %>%
  dplyr::filter(., include) %>%
  tidyr::separate(., v4_coord, c("chr", "start", "end"), remove = F) %>%
  dplyr::mutate(start = as.numeric(start) - 1, # put into bed style 0-based coordinates
                end = as.numeric(end),
                type = "gene")


# isolate just gene coordinates on v4
genes_coord <- dplyr::select(genes, chr, start, end, name_short) %>%
  dplyr::arrange(as.integer(chr), start)


# add 20kb to both ends of each gene
genes_coord_20kb = bedr(
  input = list(i = genes_coord),
  method = "slop",
  check.chr = F,
  params = paste("-g", genome_file, "-b 20000")
)

# compare to outliers in sympatric mexicana and maize
mex_maize <- c("mexicana", "maize")
genes_outliers <- list(mexicana = NULL, maize = NULL)
sig_cutoffs <- list(mexicana = NULL, maize = NULL)
for (zea in mex_maize){
  # load data
  load(paste0("local_ancestry/results/ancestry_hmm/", 
       prefix, "/Ne", Ne, "_", yesno, "Boot/anc/", 
       zea, ".pop.meta.RData")) # pop meta data
  load(paste0("local_ancestry/results/ancestry_hmm/", 
              prefix, "/Ne", Ne, "_", yesno, "Boot/anc/", 
              zea, ".pops.anc.RData")) # ancestry mean
  load(paste0("ZAnc/results/",
              prefix, "/Ne", Ne, "_", yesno, "Boot/",
              zea, ".lmElev.fit.RData")) # lmElev fits
  load(paste0("ZAnc/results/",
              prefix, "/Ne", Ne, "_", yesno, "Boot/",
              zea, ".meanAnc.fdr.RData")) # fdr ancestry mean
  assign("FDRs_meanAnc", FDRs) # give unique name
  load(paste0("ZAnc/results/",
              prefix, "/Ne", Ne, "_", yesno, "Boot/",
              zea, ".lmElev.fdr.RData")) # fdr lmElev
  assign("FDRs_lmElev", FDRs) # give unique name
  rm(FDRs) # remove ambiguous name
  FDRs <- bind_rows(mutate(FDRs_meanAnc, stat = "meanAnc"),
                    mutate(FDRs_lmElev, stat = "lmElev"))
  
  # combine lmElev fits and bed for sites
  bed_fits <- bind_cols(bed_sites, fits)
  
  # find all windows that overlap a gene of interest
  # slopes with elevation
  map_lmElev = bedr(
    input = list(a = dplyr::select(bed_fits, 
                                   c("chr", "start", "end", "(Intercept)", "envWeights")), 
                 b = genes_coord_20kb), 
    method = "map", 
    check.chr = F,
    params = paste("-g", genome_file, "-sorted -header -c 4 -o distinct")
  ) %>%
    data.table::setnames(c("chr", "start", "end", "intercept", "slope", "name_short")) %>%
    dplyr::filter(name_short != ".") %>%
    dplyr::mutate(start = as.numeric(start),
                  end = as.numeric(end),
                  intercept = as.numeric(intercept),
                  slope = as.numeric(slope))
  
  summary_lmElev <- map_lmElev %>%
    group_by(chr, name_short) %>%
    dplyr::summarise(start = min(start),
                     end = max(end),
                     min_slope = min(slope),
                     max_slope = max(slope)) %>%
    dplyr::mutate(min_sig = min_slope <= FDRs$thesholds[FDRs$tail == "low" & 
                                                          FDRs$FDR == 0.05 &
                                                          FDRs$stat == "lmElev"],
                  max_sig = max_slope >= FDRs$thesholds[FDRs$tail == "high" & 
                                                          FDRs$FDR == 0.05 &
                                                          FDRs$stat == "lmElev"]) %>%
    dplyr::mutate(min_slope = paste0(round(min_slope, 5), 
                                     ifelse(min_sig, "*", "")),
                  max_slope = paste0(round(max_slope, 5), 
                                     ifelse(max_sig, "*", "")))
  # View(summary_lmElev)
  
  # mean ancestry
  map_meanAnc = bedr(
    input = list(a = dplyr::mutate(bed_sites[ , c("chr", "start", "end")],
                                   meanAnc = anc_mean),
                 b = genes_coord_20kb), 
    method = "map", 
    check.chr = F,
    params = paste("-g", genome_file, "-sorted -header -c 4 -o distinct")
  ) %>%
    data.table::setnames(c("chr", "start", "end", "meanAnc", "name_short")) %>%
    dplyr::filter(name_short != ".") %>%
    dplyr::mutate(start = as.numeric(start),
                  end = as.numeric(end),
                  meanAnc = as.numeric(meanAnc))
  
  # what are cutoffs for low ancestry @ 2-5% of the genome?
  bed_anc <- dplyr::mutate(bed_sites, meanAnc = anc_mean) %>%
    arrange(meanAnc) %>%
    dplyr::mutate(length_Mb = length/10^6,
                  cum_length = cumsum(length_Mb),
                  percentile = cum_length/max(cum_length)*100)
  # 2-5% based on genomic length (bp) covered
  percentiles_anc <- bed_anc %>%
    summarise(perc_2 = min(meanAnc[percentile >= 2]),
              perc_5 = min(meanAnc[percentile >= 5]),
              perc_98 = max(meanAnc[percentile <= 98]),
              perc_95 = max(meanAnc[percentile <= 95]))
  
  # 2-5% of windows/tracts
  # percentiles_tracts <- bed_anc %>% 
  #   mutate(n = 1:nrow(.), percentile = n/max(n)*100) %>% 
  #   summarise(perc_2 = min(meanAnc[percentile >= 2]),
  #           perc_5 = min(meanAnc[percentile >= 5]),
  #           perc_98 = max(meanAnc[percentile <= 98]),
  #           perc_95 = max(meanAnc[percentile <= 95]))
  
  
  # add in new columns, min, max, FDR
  summary_meanAnc <- map_meanAnc %>%
    group_by(chr, name_short) %>%
    summarise(start = min(start),
              end = max(end),
              min_ancestry = min(meanAnc),
              max_ancestry = max(meanAnc)) %>%
    dplyr::mutate(min_sig = min_ancestry <= FDRs$thesholds[FDRs$tail == "low" & 
                                                             FDRs$FDR == 0.05 &
                                                             FDRs$stat == "meanAnc"],
                  max_sig = max_ancestry >= FDRs$thesholds[FDRs$tail == "high" & 
                                                             FDRs$FDR == 0.05 &
                                                             FDRs$stat == "meanAnc"],
                  min_perc2 = min_ancestry <= percentiles_anc$perc_2,
                  min_perc5 = min_ancestry <= percentiles_anc$perc_5,
                  max_perc2 = max_ancestry >= percentiles_anc$perc_98,
                  max_perc5 = max_ancestry >= percentiles_anc$perc_95) %>%
    dplyr::mutate(min_ancestry = paste0(round(min_ancestry, 5), 
                                        c("", "+", "++")[1 + min_perc2 + min_perc5], 
                                        ifelse(min_sig, "*", "")),
                  max_ancestry = paste0(round(max_ancestry, 5), 
                                        c("", "+", "++")[1 + max_perc2 + max_perc5], 
                                        ifelse(max_sig, "*", "")))
  # View(summary_meanAnc)
  
  # make table: each gene, zea, prediction, high mean anc, low mean anc, high slope, low slope
  # stars for 2%, 5% cutoffs and 5% FDR
  genes_outliers[[zea]] <- left_join(genes, summary_meanAnc, by = "name_short") %>%
    left_join(., summary_lmElev, by = "name_short") %>%
    dplyr::select(name_short, gene, category, phenotype, citation, ensembl.org_id, chr, start, end, other_notes, max_slope, min_slope, max_ancestry, min_ancestry, admixture_prediction, phenotype, v4_coord) %>%
    dplyr::mutate(subspecies = zea)
  
  # also make table with 2% and 5% cutoffs, and all FDR cutoffs
  sig_cutoffs[[zea]] <- data.frame(percentile = c(0.02, 0.05, 0.98, 0.95),
                                   threshold = c(percentiles_anc[["perc_2"]],
                                                 percentiles_anc[["perc_5"]],
                                                 percentiles_anc[["perc_98"]],
                                                 percentiles_anc[["perc_95"]]),
                                   stat = "meanAnc",
                                   tail = c("low", "low", "high", "high")) %>%
    bind_rows(., rename(FDRs, threshold = thesholds)) %>%
    arrange(stat, tail, FDR) %>%
    dplyr::select(stat, tail, FDR, percentile, threshold, n_SNPs, prop_SNPs) %>%
    mutate(subspecies = zea)
}

# combined into 1 dataframe
outliers <- left_join(dplyr::select(genes_outliers[["maize"]], -subspecies), 
                      dplyr::select(genes_outliers[["mexicana"]], name_short, max_slope, min_slope, max_ancestry, min_ancestry),
                      by = "name_short",
                      suffix = c(".maize", ".mexicana")) %>%
  dplyr::select(., c("name_short", "gene", "category", "phenotype", "admixture_prediction", "v4_coord", 
            "max_slope.maize", "min_slope.maize", "max_ancestry.maize", "min_ancestry.maize",
            "max_slope.mexicana", "min_slope.mexicana", "max_ancestry.mexicana", "min_ancestry.mexicana",
            "ensembl.org_id", "chr", "start", "end", "other_notes", "citation")) %>%
  arrange(category, chr, start)

# print raw results
write.table(outliers, file = txt_out, col.names = T, row.names = F, quote = F, sep = "\t")

# print pretty latex table
# gene, short, v4_coord, max slope\nin maize, max slope\n in mexicana, 
# min introgression\n in maize, min introgression\n in mexicana, phenotype, citation

tbl_outliers <- outliers %>%
  dplyr::select(category, name_short, v4_coord, max_slope.maize, min_ancestry.maize,
                max_slope.mexicana, max_ancestry.mexicana) %>%
  dplyr::rename(name = name_short,
                `v4 coordinates` = v4_coord,
                `max slope\nin maize` = max_slope.maize,
                `max slope\nin mexicana` = max_slope.mexicana,
                `min introgression\n in maize` = min_ancestry.maize,
                `min introgression\n in mexicana` = max_ancestry.mexicana)

print(xtable(tbl_outliers,
             type = "latex",
             latex.environments = NULL),
      include.rownames = F,
      file = tbl_out)

# print raw results
write.table(outliers, file = txt_out, col.names = T, row.names = F, quote = F, sep = "\t")

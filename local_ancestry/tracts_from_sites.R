#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)

# this script creates a bed file (set of genomic tracts)
# from a file with site position on 1 or more of the 10 maize chromosomes
# each site become one tract, and adjacent sites have tract ends that are approximately halfway between (in cM)
# the tracts together span contiguously from the first to last site on each chr
# map positions (cM) are interpolated using linear approximation and the Ogut 0.2cM genetic map
# for points outside of the original linkage map (at ends of the chr),
# we use the recombination rate from the nearest mapped window ('EXTENDED' map)

# print working directory
# setwd("~/Documents/gitErin/hilo")
print(paste("R working directory:", getwd()))

# load variables from Snakefile
rmap_file = snakemake@input[["rmap"]] # recombination map
# rmap_file = "data/linkage_map/ogut_2015_rmap_v2_to_v4_EXTENDED.txt"
sites_file = snakemake@input[["sites"]]
# sites_file = "local_ancestry/results/thinnedSNPs/HILO_MAIZE55/K2/whole_genome.var.sites"
bed_file = snakemake@output[["bed"]]
# bed_file = "local_ancestry/results/thinnedSNPs/HILO_MAIZE55/K2/whole_genome.bed"

# load sites file
sites = read.table(sites_file, header = F, sep = "\t") %>%
  data.table::setnames(., c("chr", "pos", "allele1", "allele2"))
# split into a list by chromosome
sites_by_chr = lapply(1:10, function(i) dplyr::filter(sites, chr == i))

# load recombination map for the whole genome
rmap_genome = read.table(rmap_file, header = T, sep = "\t", stringsAsFactors = F)
# divide into 10 maps -- 1 per chromosome
rmap_by_chr = lapply(1:10, function(i) dplyr::filter(rmap_genome, chr == i))

# find tract start and end bp positions,
# going chromosome by chromosome
for (i in 1:10){
  if (nrow(sites_by_chr[[i]]) > 0){ # if there are any sites on that chr
    # get cM position for each site by linear interpolation
    pos_cM = stats::approx(x = rmap_by_chr[[i]]$pos_bp, y = rmap_by_chr[[i]]$pos_cM,
                                     xout = sites_by_chr[[i]]$pos, method = "linear")$y
    # get midpoints between sites (starts) in cM for each tract
    start_cM = pos_cM[1:(length(pos_cM) - 1)] + diff(pos_cM, lag = 1)/2
    # find nearest bp (round to integer) for each cM start position
    start_bp = round(stats::approx(x = rmap_by_chr[[i]]$pos_cM, y = rmap_by_chr[[i]]$pos_bp,
                            xout = start_cM, method = "linear")$y)
    sites_by_chr[[i]]$start = c(sites_by_chr[[i]]$pos[1], start_bp) - 1 # minus 1 because bed files are zero indexed
    sites_by_chr[[i]]$end = c(start_bp - 1, sites_by_chr[[i]]$pos[length(pos_cM)])
  }
}

# put all chromosomes back together into 1 bed file
bed <- do.call(rbind, sites_by_chr) %>%
  dplyr::select(chr, start, end, pos)

# write output bed file
write.table(bed, bed_file,
            col.names = F, row.names = F,
            sep = "\t", quote = F)

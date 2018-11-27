#!/usr/bin/env Rscript

# this script extends rmapALL genetic map to the first and last bp
# of each chromosome, using furthest left/right map rate
source("calcMapPosCM.R")

# load original recombination map
# cleaned for positions out of order 
# after lifting onto genome version 4 ("INCLUDED")
rmapALL <- read.table("../data/linkage_map/ogut_fifthcM_map_agpv4_INCLUDE.txt",
                      sep = "\t", header = F, stringsAsFactors = F)
colnames(rmapALL) <- c("SNP", "marker", "pos_cM", "chr", "pos_bp")

# load chromosome ranges
zea_chr <- read.table("../data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths",
                      sep = "\t", header = F, stringsAsFactors = F)
colnames(zea_chr) <- c("chr", "length")

# last position each chromosome:
zea_chr$map_ends = apply(zea_chr, 1, function(i) 
  calcMapPosCM(chr = i["chr"], pos_bp = i["length"],
               rmap = rmapALL))
# first position each chromosome:
zea_chr$map_starts = apply(zea_chr, 1, function(i)
  calcMapPosCM(chr = i["chr"], pos_bp = 1,
               rmap = rmapALL))
# add start and end of chr genetic map positions to rmap
rmapEXT_unordered = rbind(rmapALL,
                          data.frame(SNP = paste0("start_chr", 
                                                  zea_chr$chr), 
                                     marker = paste0("Start", zea_chr$chr),
                                     chr = zea_chr$chr, 
                                     pos_cM = zea_chr$map_starts,
                                     pos_bp = 1),
                          data.frame(SNP = paste0("end_chr", 
                                                  zea_chr$chr), 
                                     marker = paste0("End", zea_chr$chr),
                                     chr = zea_chr$chr, 
                                     pos_cM = zea_chr$map_ends,
                                     pos_bp = zea_chr$length))
# order by chromosome, then by bp
rmapEXT = rmapEXT_unordered[
  with(rmapEXT_unordered, order(chr, pos_bp)),
  ]
write.table(rmapEXT, "../data/linkage_map/ogut_fifthcM_map_agpv4_EXTENDED.txt",
            sep = "\t", col.names = T, row.names = F, quote = F)

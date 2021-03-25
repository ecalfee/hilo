#!/usr/bin/env Rscript
# Author: Erin Calfee. Last updated 03/2021.

# this script extends rmap_included genetic map to the first and last bp
# of each chromosome, using furthest left/right map rate

# load input and output file names with snakemake 
# (or alternatively use commented out lines to run outside of a snakemake pipeline)
rmap_clean = snakemake@input[["rmap_clean"]]
# rmap_clean = "linkage_map/results/ogut_2015_rmap_v2_to_v4_INCLUDED.txt"
genome = snakemake@input[["genome"]]
# genome = "data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"
rmap_ext = snakemake@output[["rmap_ext"]]
# rmap_ext = "linkage_map/results/ogut_2015_rmap_v2_to_v4_EXTENDED.txt"

# load original recombination map
# cleaned for positions out of order 
# after lifting onto genome version 4 ("INCLUDED")
rmap_included <- read.table(rmap_clean,
                      sep = "\t", header = T, stringsAsFactors = F)

# load chromosome ranges
zea_chr <- read.table(genome,
                      sep = "\t", header = F, stringsAsFactors = F)
colnames(zea_chr) <- c("chr", "length") # chr = chromosomes 1-10, length = chromosome lengths

# function to calculate map position in cM of any SNP in bp
calcMapPosCM = function(chrom, pos_bp, rmap){
  # subset recombination map to only relevant chromosome
  rmapChr = rmap[rmap$chr == chrom, ]
  # if position is outside the range of the map 
  # (less than 1st value or greater than last value)
  # use avg. recombination rate for the closest bin (at that end of the chr)
  if (pos_bp < rmapChr$pos_bp[1]){
    #print("off map left")
    left = rmapChr[1, ] # leftmost mapped position
    right = rmapChr[2, ] # next position from furthest left
  } else{
    if (pos_bp >= rmapChr$pos_bp[nrow(rmapChr)]){
      left = rmapChr[nrow(rmapChr)-1, ] # second to last position from the right
      right = rmapChr[nrow(rmapChr), ] # rightmost mapped position on chromosome
    } else{ # otherwise if SNP is within the map, find the closest mapped positions above and below pos
      left = rmapChr[rmapChr["pos_bp"] <= pos_bp, ][sum(rmapChr$pos_bp <= pos_bp), ] # last mapped position smaller than or equal to pos
      right = rmapChr[rmapChr["pos_bp"] > pos_bp, ][1, ] # first mapped position greater than pos
    }
  }
  # calculate cM/bp in local region between 2 closest mapped positions
  rate = (right$pos_cM - left$pos_cM) / (right$pos_bp - left$pos_bp)
  pos_cM = left$pos_cM + rate*(pos_bp - left$pos_bp)
  return(pos_cM)
}


# last position each chromosome:
zea_chr$map_ends = apply(zea_chr, 1, function(i) 
  calcMapPosCM(chr = i["chr"], pos_bp = i["length"],
               rmap = rmap_included))

# first position each chromosome:
zea_chr$map_starts = apply(zea_chr, 1, function(i)
  calcMapPosCM(chr = i["chr"], pos_bp = 1,
               rmap = rmap_included))

# add start and end of chr genetic map positions to rmap
rmap_extended_unordered = bind_rows(rmap_included,
                          data.frame(marker = paste0("Start", zea_chr$chr),
                                     chr = zea_chr$chr, 
                                     pos_cM = zea_chr$map_starts,
                                     pos_bp = 1),
                          data.frame(marker = paste0("End", zea_chr$chr),
                                     chr = zea_chr$chr, 
                                     pos_cM = zea_chr$map_ends,
                                     pos_bp = zea_chr$length)) %>%
  mutate(SNP = paste(chr, pos_bp, sep = "_")) %>%
  dplyr::select(SNP, marker, chr, pos_cM, pos_bp)

# order by chromosome, then by bp
rmap_extended = rmap_extended_unordered %>%
  arrange(chr, pos_bp)

# write output to file
write.table(rmap_extended, file = rmap_ext,
            sep = "\t", col.names = T, row.names = F, quote = F)

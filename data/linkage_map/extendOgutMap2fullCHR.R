#!/usr/bin/env Rscript

# this script extends rmapALL genetic map to the first and last bp
# of each chromosome, using furthest left/right map rate

# load original recombination map
# cleaned for positions out of order 
# after lifting onto genome version 4 ("INCLUDED")
rmapALL <- read.table("ogut_fifthcM_map_agpv4_INCLUDE.txt",
                      sep = "\t", header = F, stringsAsFactors = F)
colnames(rmapALL) <- c("SNP", "marker", "pos_cM", "chr", "pos_bp")

# load chromosome ranges
zea_chr <- read.table("../data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths",
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
rmapEXT = rmapEXT_unordered[with(rmapEXT_unordered, order(chr, pos_bp)), ]
write.table(rmapEXT, "ogut_fifthcM_map_agpv4_EXTENDED.txt",
            sep = "\t", col.names = T, row.names = F, quote = F)

# useful functions for recombination map

# load recombination map
rmapALL <- read.table("../data/linkage_map/ogut_fifthcM_map_agpv4_INCLUDE.txt",
                  sep = "\t", header = F, stringsAsFactors = F)
colnames(rmapALL) <- c("SNP", "marker", "pos_cM", "chr", "pos_bp")

# load chromosome ranges
zea_chr <- read.table("../data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths",
                     sep = "\t", header = F, stringsAsFactors = F)
colnames(zea_chr) <- c("chr", "length")

# function to calculate map rate of any SNP (also implemented in calcMapPos.py)
calcMapRate = function(chrom, pos, rmap = rmapALL){
  # subset recombination map to only relevant chromosome
  rmapChr = rmap[rmap$chr == chrom, ]
  # if position is outside the range of the map 
  # (less than 1st value or greater than last value)
  # use avg. recombination rate for the closest bin (at that end of the chr)
  if (pos < rmapChr$pos_bp[1]){
    #print("off map left")
    left = rmapChr[1, ] # leftmost mapped position
    right = rmapChr[2, ] # next position from furthest left
  } else{
    if (pos >= rmapChr$pos_bp[nrow(rmapChr)]){
      left = rmapChr[nrow(rmapChr)-1, ] # second to last position from the right
      right = rmapChr[nrow(rmapChr), ] # rightmost mapped position on chromosome
    } else{ # otherwise if SNP is within the map, find the closest mapped positions above and below pos
      left = rmapChr[rmapChr["pos_bp"] <= pos, ][sum(rmapChr$pos_bp <= pos), ] # last mapped position smaller than or equal to pos
      right = rmapChr[rmapChr["pos_bp"] > pos, ][1, ] # first mapped position greater than pos
    }
  }
  # calculate cM/bp in local region between 2 closest mapped positions
  rate = (right$pos_cM - left$pos_cM) / (right$pos_bp - left$pos_bp)
  return(rate)
}

# function to calculate map position in cM of any SNP in bp
# uses map positions in rmap as a baseline (can be negative!)
calcMapPosCM = function(chrom, pos_bp, rmap = rmapALL){
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

# function to calculate bp position for an cM position input
# uses map positions in rmap as a baseline (can be negative!)
# note if round = F, pos_bp could be a fraction -- not rounded
calcMapPosBP = function(chrom, pos_cM, rmap = rmapALL, zea = zea_chr, round = T){
  # subset recombination map to only relevant chromosome
  rmapChr = rmap[rmap$chr == chrom, ]
  # if position is outside the range of the map 
  # (less than 1st value or greater than last value)
  # use avg. recombination rate for the closest bin (at that end of the chr)
  if (pos_cM < rmapChr$pos_cM[1]){
    #print("off map left")
    left = rmapChr[1, ] # leftmost mapped position
    right = rmapChr[2, ] # next position from furthest left
  } else{
    if (pos_cM >= rmapChr$pos_cM[nrow(rmapChr)]){
      left = rmapChr[nrow(rmapChr)-1, ] # second to last position from the right
      right = rmapChr[nrow(rmapChr), ] # rightmost mapped position on chromosome
    } else{ # otherwise if SNP is within the map, find the closest mapped positions above and below pos
      left = rmapChr[rmapChr["pos_cM"] <= pos_cM, ][sum(rmapChr$pos_cM <= pos_cM), ] # last mapped position smaller than or equal to pos
      right = rmapChr[rmapChr["pos_cM"] > pos_cM, ][1, ] # first mapped position greater than pos
    }
  }
  # calculate cM/bp in local region between 2 closest mapped positions
  rate = (right$pos_cM - left$pos_cM) / (right$pos_bp - left$pos_bp)
  pos_bp = left$pos_bp + (1/rate)*(pos_cM - left$pos_cM)
  if (round) pos_bp = round(pos_bp) # round to nearest whole number bp position
  
  # correct found special boundary cases:
  if (pos_bp < 0){ # off chromosome to the left
    pos_bp = 0 
    pos_cM = calcMapPosCM(chrom = chrom, pos_bp = pos_bp, rmap = rmapChr)
  }
  end_chr = zea[zea$chr==chrom, "length"] # end of chromosome
  if (pos_bp > end_chr){ # off chromosome to the right
    pos_bp = end_chr
    pos_cM = calcMapPosCM(chrom = chrom, pos_bp = pos_bp, rmap = rmapChr)
  }
  
  return(c(pos_bp = pos_bp, pos_cM = pos_cM, rate = rate))
  # needs to deal with being negative positions bp --> 0 
  # or exceeding length of chromosome (needs chromosome lengths as input)
}


# test functions for a few positions
#calcMapRate(chr = 1 , pos = 412)
#calcMapRate(chr = 1 , pos = 1000)
#calcMapPosCM(chr = 2, pos_bp = 500000, rmap = rmapALL)
#calcMapPosCM(chr = 2, pos_bp = 5000000)
#calcMapPosBP(chr = 2, pos_cM = 2, round = T)
#calcMapPosBP(chr = 2, pos_cM = 109, round = F)
#calcMapPosBP(chr = 2, pos_cM = 8888888, round = F)
#calcMapPosBP(chr = 2, pos_cM = -70, round = F)


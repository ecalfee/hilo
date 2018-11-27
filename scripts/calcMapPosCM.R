# function to calculate map position in cM of any SNP in bp
# uses map positions in rmap as a baseline (can be negative!)
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

# load recombination map and make a function to calculate map position of any SNP
rmap = read.table("../data/linkage_map/ogut_fifthcM_map_agpv4_INCLUDE.txt",
                  sep = "\t", header = F, stringsAsFactors = F)
colnames(rmap) = c("SNP", "marker", "pos_cM", "chr", "pos_bp")

# helper function to calculate map position (also implemented in calcMapPos.py)
calcMapRate = function(chrom, pos, rmap){
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
dir_in = "../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/input/"
pos_thin <- read.table(paste0(dir_in, "pop18.anc_hmm.input"),
                       stringsAsFactors = F)[ , 1:2]
colnames(pos_thin) <- c("chr", "pos")
# test function for a few positions
#calcMapRate(chr = 1 , pos =412, rmap = rmap)
#calcMapRate(chr = 1 , pos =1000, rmap = rmap)
# calculate map rate for every position in pos_thin
pos_thin$rate <- apply(pos_thin, 1, 
                       function(i) calcMapRate(chr = i["chr"], 
                                               pos = i["pos"],
                             rmap = rmap))
write.table(pos_thin, 
            paste0(dir_in, "pos_recomb_rates.txt"), 
            col.names = T, row.names = F, quote = F, sep = "\t")


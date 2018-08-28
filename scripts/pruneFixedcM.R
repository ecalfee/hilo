# THIS IS A DRAFT SCRIPT -- NOT NEEDED AT THE MOMENT BECAUSE OF 
# pruneFixedcM.py
# inputs (e.g. from last processed file & other variables for thinning):
input_kept_cM = 0
input_kept_chr = 0 # not a true chromosome
minDist = .0001 #(.001cM ~ 1Kb on average across the genome)

regions = region_start:region_end # takes in consecutive set of regions to apply thinning
region_start = 9
region_end = 9

# load recombination map and make a function to calculate map position of any SNP
rmap = read.table("../data/linkage_map/ogut_fifthcM_map_agpv4_INCLUDE.txt",
                  sep = "\t", header = F, stringsAsFactors = F)
colnames(rmap) = c("SNP", "marker", "pos_cM", "chr", "pos_bp")

# helper function to calculate map position (also implemented in calcMapPos.py)
calcMapPos = function(chrom, pos, rmap){
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
  #print("cM/bp rate is " + str(rate)) 
  # use this rate to find distance from leftmost position
  # note that the bp difference * rate will be negative if the position is off the map to the left of all mapped positions
  map_pos = left$pos_cM + (pos - left$pos_bp)*rate 
  
  if (rate <= 0){ # should always be true with cleaned recombination rate map
    stop(paste("Recomb. rate calculated is NOT strictly positive, pos = ", pos))
  } 
  return(map_pos) # returns position in cM (note: may be negative! b/c mapped positions go negative from some previous 0)
}

# function to thin SNPs based on a minimum cM distance between SNPs
thin_SNPs = function(d, rmap, D, minInd, minDist, first_cM, first_chr){
  last_kept_cM = first_cM
  last_kept_chr = first_chr
  # now filter SNPs for LD
  # saves last kept SNP and a vector of T/F for which SNPs are to be kept
  keep = rep(F, nrow(d))
  d1 = d %>%
    mutate(., cM = sapply(1:nrow(.), function(i) 
      calcMapPos(chrom = chromo[i], pos = position[i], rmap = rmap)))
  for (i in 1:length(keep)){
    # always in low LD if on diff. chromosomes
    cM1 = d1[i, "cM"]
    chr1 = d1[i, "chromo"]
    keep[i] = abs(cM1 - last_kept_cM) >= minDist | chr1 != last_kept_chr
    if (keep[i]){ # if you keep SNP, update last-kept-SNP = current NSP
      last_kept_cM = cM1
      last_kept_chr = chr1
    }
  }
  d_thin = d1[keep, ]
  # calculate diference in morgans between positions kept
  d_thin$diff_M = diff(c(first_cM, d_thin$cM))/100
  ifelse(first_chr != last_kept_chr){ # if this file is the first of a new chromosome
    d_thin$diff_M[1] = 1 # make Morgan difference = 1 (arbitrary for ancestry_hmm input)
  }
  return(d_thin)
}


# print out a file with chromosome, position, major, minor & difference between adjacent 
# SNPs in Morgans (to 6 or 8 decimal places .. I'd have to check)
# I need to ultimately get 1 file per chromosome after thinning .. but I could concatenate after thinning
# or I could still concatenate at the very end after getting frequencies etc.
last_cM = 0
last_chr = 0
for (r in regions){
  
  d = thin_SNPs(d = getSNPs(region = r, 
                            dir_gl = dir_gl, 
                            dir_sites = dir_sites, 
                            allo_maize = allo_maize, 
                            allo_mex = allo_mex), 
                rmap = rmap, 
                D = D, 
                minInd = minInd, 
                minDist = minDist, 
                first_cM = last_cM, 
                first_chr = last_chr)
  # write to file
  write.table(paste0(dir_sites, "/thinned/")) ######### NOT FINISHED
}



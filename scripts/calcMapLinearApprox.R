# this script creates the functions to calculate the
# approximate cM or bp position for any site on a chromosome
# using linear interpolation
# e.g. bp2cM(chrom = 1, pos_bp = 1) and cM2bp(chrom = 5, pos_cM = 0.2)
# These functions return NA if the position is out of bounds (as in, off the chromosome)

# load extended recombination map (extended to full length of chromosomes)
rmapEXT = read.table("../data/linkage_map/ogut_2015_rmap_v2_to_v4_EXTENDED.txt",
                     sep = "\t", stringsAsFactors = F, header = T)


# testing approx() linear interpolation function
# returns NA if outside physical bp range of chromosome
bp2cM_byChr = function(chrom, rmap = rmapEXT){
  rmapCHR = rmap[rmap$chr == chrom,]
  approxfun(x = rmapCHR$pos_bp, y = rmapCHR$pos_cM)
}
cM2bp_byChr = function(chrom, rmap = rmapEXT){
  rmapCHR = rmap[rmap$chr == chrom,]
  approxfun(x = rmapCHR$pos_cM, y = rmapCHR$pos_bp)
}
# make a list of functions (1 per chromosome)
bp2cM_ALLChr = lapply(1:10, function(i)
  bp2cM_byChr(chrom = i, rmap = rmapEXT))
cM2bp_ALLChr = lapply(1:10, function(i)
  cM2bp_byChr(chrom = i, rmap = rmapEXT))

bp2cM = function(chrom, pos_bp, list_FUN = bp2cM_ALLChr){
  return(list_FUN[[chrom]](pos_bp))
}

cM2bp = function(chrom, pos_cM, list_FUN = cM2bp_ALLChr){
  return(round(list_FUN[[chrom]](pos_cM))) # round to nearest bp
}

# test cases:
#bp2cM(chrom = 1, pos_bp = 1) # -5.103104
#bp2cM(chrom = 1, pos_bp = 30) # slightly closer to 0
#bp2cM(chrom = 0, pos_bp = 1) # error
#bp2cM(chrom = 10, pos_bp = 160000000) # NA -- out of bounds
#bp2cM(chrom = 10, pos_bp = -3) # NA -- out of bounds

#cM2bp(chrom = 1, pos_cM = bp2cM(chrom = 1, pos_bp = 1)) # 1
#cM2bp(chrom = 7, pos_cM = 0) # some number
#cM2bp(chrom = 1, pos_cM = -6) # NA -- out of bounds
#cM2bp(chrom = 4, pos_cM = 300) # NA -- out of bounds

# testing subfunctions:
# pb2cM_chr1 = bp2cM_byChr(chrom = 1)
# pb2cM_chr1(c(0, 1, 50000, 300000, 400000,
#             100000000, 8000000000)) # test cases
# bp2cM_ALLChr[[1]](c(0, 1, 50000, 300000, 400000, 100000000, 8000000000))
# cM2bp_chr8 = cM2bp_byChr(chrom = 8)
# cM2bp_chr8(c(-90, -3.1, -1, 0.01, 4, 33))
# cM2bp_ALLChr[[8]](c(-90, -3.1, -1, 0.01, 4, 33))
# cM2bp_chr8(seq(-1, 1, by = .01)) # does produce decimals
# round output to nearest bp
# round(cM2bp_ALLChr[[8]](seq(-1, 1, by = .01)))

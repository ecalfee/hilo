library(IRanges)
source("rmap_functions.R")

# this file has two functions for gene density:
# (1) densityCDS = number of coding BP/ total # bp in a fixed cM window
# around a site
# (2) linkageCDS = summing over any contigous CDS
# that overlap at all with a fixed window around a site
# what is the cumulative # bp CDS * cM distance away

# CoDing Sequences: CDS have been merged into contiguous ranges
cds_genome = lapply(1:10, function(i) 
read.table(paste0("../data/refMaize/geneAnnotations/CDS_", i, ".txt"),
           sep = "\t", header = T, stringsAsFactors = F))

# local density of CoDing Sequences / bp
# wind = window size in cM
densityCDS = function(chrom, pos, wind = 1, cds = cds_genome){
  # position focal locus
  pos_cM = calcMapPosCM(chrom = chrom, pos_bp = pos)
  # left and right boundaries of window
  left = calcMapPosBP(chrom = chrom, pos_cM = pos_cM - .5*wind)
  right = calcMapPosBP(chrom = chrom, pos_cM = pos_cM + .5*wind)
  
  # number of bp overlapping with CDS ranges
  ir = IRanges(start = cds[[chrom]]$start,
               end = cds[[chrom]]$end)
  my_wind = IRanges(start = left["pos_bp"], end = right["pos_bp"])
  hits = findOverlaps(my_wind, ir) # which CDS are within my window (at least partially)?
  overlap = pintersect(my_wind[queryHits(hits)], 
             ir[subjectHits(hits)]) # what is the length of overlap?
  bpCDS = sum(width(overlap)) # sum all overlap, i.e. CDS bp within window
  bpTOT = width(my_wind)
  return(c(dens = bpCDS/bpTOT, bpCDS = bpCDS, bpTOT = bpTOT, widthCM = as.numeric(right["pos_cM"] - left["pos_cM"])))
}
# test function
#densityCDS(chrom = 1, pos = 0, wind = .1) # edge
#densityCDS(chrom = 1, pos = 10000, wind = 1) # still on edge
#densityCDS(chrom = 1, pos = 600000, wind = 1) # full cM, no edge effects

# score for linkage to CoDing Sequences
# default = 0 for no CDS in a window
linkageCDS = function(chrom, pos, wind = 1, cds = cds_genome){

  # position focal locus
  pos_cM = calcMapPosCM(chrom = chrom, pos_bp = pos)
  # left and right boundaries of window
  left = calcMapPosBP(chrom = chrom, pos_cM = pos_cM - .5*wind)
  right = calcMapPosBP(chrom = chrom, pos_cM = pos_cM + .5*wind)

  # overlapping CDS ranges
  ir = IRanges(start = cds[[chrom]]$start,
               end = cds[[chrom]]$end)
  my_wind = IRanges(start = left["pos_bp"], end = right["pos_bp"])
  widthCM = as.numeric(right["pos_cM"] - left["pos_cM"])
  bpTOT = width(my_wind)
  
  # which CDS are within my window (at least partially)?
  hits = ir[subjectHits(findOverlaps(my_wind, ir))] 
  
  if (length(hits) == 0){ # no CDS overlap with window
    return(c(linkageCDS = 0, bpTOT = bpTOT, widthCM = widthCM))
  } else{
    # find midpoint of each overlapping CDS
    midpoints = (start(hits) + end(hits))/2
    # get distance, in cM, from focal site for each CDS midpoint
    cM_dist = sapply(midpoints, function(x) 
      abs(calcMapPosCM(chrom = chrom, pos_bp = x) - pos_cM))
    # sum the genetic distance from the focal locus
    # times the # bp in each overlapping CoDing Sequence
    linkage = sum(cM_dist*width(hits))
    return(c(linkageCDS = linkage, bpTOT = bpTOT,
             widthCM = widthCM))
  }
}
# test function
linkageCDS(chrom = 1, pos = 0, wind = .001) # edge
linkageCDS(chrom = 1, pos = 500000, wind = .01) # edge
linkageCDS(chrom = 1, pos = 5000000, wind = 1)

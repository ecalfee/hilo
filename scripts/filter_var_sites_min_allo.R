library(dplyr)
# this script filters var sites before an LD filter
# for ultimately inputting into ancestry_hmm for local ancestry calling

minInd = 4 # minimum number of individuals per population to include a site
D = 0.3 # minimum difference in estimated allele frequency between 
# allopatric maize and allopatric mexicana (goal = to enrich for informative sites)

# prefix/main directory for data
dir_gl = "../data/geno_lik/merged_pass1_all_alloMaize4Low_16/allVar/"
dir_sites = "../data/var_sites/merged_pass1_all_alloMaize4Low_16/"
region = "region_9"
allo_maize = "maize.allo.4Low16"
allo_mex = "mexicana.allo.withXochi35"

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


# get allopatric MAF data for some region
maize = read.table(gzfile(paste0(dir_gl, allo_maize, "/", 
                                 region, ".mafs.gz")), header=T, stringsAsFactors = F) %>%
  rename(.,  MAF_maize = knownEM) %>%
  rename(., nInd_maize = nInd) %>%
  select(., -ref) # remove reference

mex = read.table(gzfile(paste0(dir_gl, allo_mex, "/",
                               region, ".mafs.gz")), header=T, stringsAsFactors = F) %>%
  rename(.,  MAF_mex = knownEM) %>%
  rename(., nInd_mex = nInd) %>%
  select(., -ref) # remove reference

#d2 = read.table(gzfile(paste0(dir_gl, "whole_genome_pruned_every_1000.beagle.gz")), header=T, stringsAsFactors = F)
#sum(duplicated(d2[ , 1:4])) # no duplicated sites by chance in the every 1000 SNPs file

d = read.table(paste0(dir_sites, region, ".var.sites"), header=F, stringsAsFactors = F) %>%
  setNames(., colnames(maize)[1:4]) %>%
  left_join(., maize, by = colnames(maize)[1:4]) %>%
  left_join(., mex, by = colnames(maize)[1:4]) %>%
  filter(., !duplicated(.)) # temporary problem with duplicated positions -- I'll fix upstream later

# thin markers for informative sites
d_thin = d %>% 
  filter(., nInd_maize >= minInd & nInd_mex >= minInd) %>%
  filter(., abs(MAF_maize - MAF_mex) >= D) %>%
  mutate(., cM = sapply(1:nrow(.), function(i) 
    calcMapPos(chrom = chromo[i], pos = position[i], rmap = rmap)))


# try different filtering and look at results
# what % of sites will be left after some minor filtering for more informative sites?
tot = sum(complete.cases(d)) # number of sites with data from maize and mex (largest possible total)
nrow(d_thin)/tot
# at a filter of at least 4 individuals per allopatric population,
# around .2 allele freq. difference or .3 is ~30% and ~17% of sites remaining post-filtering
# D = .4 drops the number of sites down to 1% 
# so my estimate is D = .3 will be best for ancestry calling to still get high enough
# density of SNPs but enriching for more informative SNPs .. I will still need to do some 
# filtering for SNPs in very high LD (too close cM distance)

# now filter SNPs for LD

# inputs (e.g. from last processed file):
input_kept_cM = 0
input_kept_chr = 0 # not a true chromosome
minDist = .0001 #(.001cM ~ 1Kb on average across the genome)

# saves last kept SNP and a vector of T/F for which are to be kept
keep = rep(F, nrow(d_thin))
last_kept_cM = input_kept_cM
last_kept_chr = input_kept_chr 

for (i in 1:length(keep)){
  # always in low LD if on diff. chromosomes
  cM1 = d_thin[i, "cM"]
  chr1 = d_thin[i, "chromo"]
  keep[i] = abs(cM1 - last_kept_cM) >= minDist | chr1 != last_kept_chr
  if (keep[i]){ # if you keep SNP, update last-kept-SNP = current NSP
    last_kept_cM = cM1
    last_kept_chr = chr1
  }
}
d_thinner = d_thin[keep,]
#sum(keep)
#hist(d_thinner$nInd_maize)
#hist(d_thinner$nInd_mex)
#table(d_thinner$nInd_mex)

# print out a file with chromosome, position, major, minor & difference between adjacent 
# SNPs in Morgans (to 6 or 8 decimal places .. I'd have to check)
# I need to ultimately get 1 file per chromosome after thinning .. but I could concatenate after thinning
# or I could still concatenate at the very end after getting frequencies etc.

library(dplyr)
# this script filters var sites before an LD filter
# for ultimately inputting into ancestry_hmm for local ancestry calling

minInd = 4 # minimum number of individuals per population to include a site
D = .3 # minimum difference in estimated allele frequency between 
# allopatric maize and allopatric mexicana (goal = to enrich for informative sites)

# load recombination map and make a function to calculate map position of any SNP
rmap = 
calcMapPos = function(chrom, pos, rmap):
  # subset recombination map to only relevant chromosome
  rmapChr = rmap[rmap["chrom"] == chrom]
# if position is outside the range of the map 
# (less than 1st value or greater than last value)
# use avg. recombination rate for the closest bin (at that end of the chr)
if (pos < rmapChr["pos_bp"].iloc[0]):
  #print("off map left")
  left = rmapChr.iloc[0, :] # leftmost mapped position
right = rmapChr.iloc[1, :] # next position from furthest left
elif (pos >= rmapChr["pos_bp"].iloc[-1]):
  left = rmapChr.iloc[-2, :] # second to last position from the right
right = rmapChr.iloc[-1, :] # rightmost mapped position on chromosome
# otherwise if SNP is within the map, find the closest mapped positions above and below pos
else:
  left = rmapChr[rmapChr["pos_bp"] <= pos].iloc[-1] # last mapped position smaller than or equal to pos
right = rmapChr[rmapChr["pos_bp"] > pos].iloc[0] # first mapped position greater than pos

# calculate cM/bp in local region between 2 closest mapped positions
rate = (right["pos_cM"] - left["pos_cM"]) / (right["pos_bp"] - left["pos_bp"])
#print("cM/bp rate is " + str(rate)) 
# use this rate to find distance from leftmost position
# note that the bp difference * rate will be negative if the position is off the map to the left of all mapped positions
map_pos = left["pos_cM"] + (pos - left["pos_bp"])*rate 

if rate <= 0: # should always be true with cleaned recombination rate map
  raise ValueError("Recomb. rate calculated is NOT strictly positive, pos = " + str(pos))

return(map_pos) # returns position in cM (note: may be negative! b/c mapped positions go negative from some previous 0)





# prefix/main directory for data
dir_gl = "../data/geno_lik/merged_pass1_all_alloMaize4Low_16/allVar/"
dir_sites = "../data/var_sites/merged_pass1_all_alloMaize4Low_16/"
region = "region_9"
allo_maize = "maize.allo.4Low16"
allo_mex = "mexicana.allo.withXochi35"

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

# try different filtering and look at results
tot = sum(complete.cases(d)) # number of sites with data from maize and mex (largest possible total)
d %>% # what % of sites will be left after some minor filtering for more informative sites?
  filter(., nInd_maize >= minInd & nInd_mex >= minInd) %>%
  filter(., abs(MAF_maize - MAF_mex) >= D) %>%
  nrow(.)/tot
# at a filter of at least 4 individuals per allopatric population,
# around .2 allele freq. difference or .3 is ~30% and ~17% of sites remaining post-filtering
# D = .4 drops the number of sites down to 1% 
# so my estimate is D = .3 will be best for ancestry calling to still get high enough
# density of SNPs but enriching for more informative SNPs .. I will still need to do some 
# filtering for SNPs in very high LD (too close cM distance)


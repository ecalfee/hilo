library(dplyr)
library(tidyr)
library(ggplot2)
# calculate diversity from the folded SFS:
popN <- c(18:31,33:35,360:363,365:374,1000,2000,3000,5000)
calcTheta <- function(popN){
  sfs <- t(read.table(paste0("../data/SFS/pass1/N1000.L100.regions/pop", popN, ".folded.sfs")))
  # SFS gives ML # of sites in each category of the folded SFS
  # actually I accidentally didn't create a folded SFS, but it's not polarized with an outgroup
  # so we can't identify the derived allele (doesn't matter)
  # 0 minor alleles, 1 minor allele, 2 minor alleles etc.
  
  sfs_perc <- sfs/sum(sfs) # percent of sites falling into each bin of SFS; divide by total number of sites
  N = length(sfs)-1
  sfs_freq <- 0:N/N # divide by total # of individuals in sample
  # heterozygosity (pi): mean across all sites of 2*p*(1-p)
  het = sum(2*sfs_freq*(1-sfs_freq)*sfs_perc)
  thetaW_num = sum(sfs[2:N]) # watterson's theta numerator = # segregating sites (ignore 1st and last bin)
  thetaW_den = sum(1/(1:(N-1))) # denominator watterson's theta is neutral branch lengths
  thetaW = thetaW_num/thetaW_den/sum(sfs) # divide by total # loci to get a per bp theta estimate
  return(c(het=het, thetaW=thetaW, N=N))
}
# calculate diversity and add population labels
# pop labels
hilo <- read.table("../data/HILO_IDs_cov_pass1.csv", stringsAsFactors = F, header = T)
groups = data.frame(popN = c(1000, 2000, 3000, 5000), zea = c("maize", "mexicana", "mexicana", "maize"),
                    symp_allo = c("sympatric", "sympatric", "allopatric", "allopatric"),
                    stringsAsFactors = F)

diversity <- data.frame(t(sapply(popN, calcTheta))) %>%
  mutate(., popN = popN) %>%
  left_join(., union(unique(hilo[ , c("popN", "zea", "symp_allo")]), groups), by = "popN")

write.table(format(diversity, digits = 4), "../data/SFS/pass1/N1000.L100.regions/diversity_pi_theta_from_SFS.txt",
            quote = F, row.names = F, sep = "\t")

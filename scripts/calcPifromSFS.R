library(dplyr)
library(tidyr)
library(ggplot2)
source("calcTheta.R") # function to calculate pi and watterson's theta diversity statistics

# calculate diversity from the folded SFS:
#popN <- c(18:31,33:35,360:363,365:374,1000,2000,3000,5000)
popN <- c(1000,2000,3000,5000)
unfolded_sfs_files <- paste0("../data/SFS/pass1/N1000.L100.regions/pop", popN, ".folded.sfs")
# despite the name, actually I accidentally didn't create a folded SFS, 
# but instead an unfolded SFS, but it's not polarized with an outgroup
# so we can't identify the derived allele (doesn't matter)
# 0 minor alleles, 1 minor allele, 2 minor alleles, N minor alleles etc.


# calculate diversity and add population labels
# pop labels
hilo <- read.table("../data/pass1_ids.txt", sep = "\t", header = T, stringsAsFactors = F)
groups = data.frame(popN = c(1000, 2000, 3000, 5000), zea = c("maize", "mexicana", "mexicana", "maize"),
                    symp_allo = c("sympatric", "sympatric", "allopatric", "allopatric"),
                    stringsAsFactors = F)

diversity <- data.frame(t(sapply(unfolded_sfs_files, calcTheta))) %>%
  mutate(., popN = popN) %>%
  left_join(., union(unique(hilo[ , c("popN", "zea", "symp_allo")]), groups), by = "popN")

write.table(format(diversity, digits = 4), "../data/SFS/pass1/N1000.L100.regions/diversity_pi_theta_from_SFS.txt",
            quote = F, row.names = F, sep = "\t")

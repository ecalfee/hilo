library(dplyr)
# this script takes in a set of sites for a region,
# and a set of populations with within-ancestry allele frequencies
# and calculates Fst for that region between all pairs

# e.g. look at chr4 inversion:
inv <- read.table("../data/refMaize/inversions/knownInv_v4_coord.txt", stringsAsFactors = F)
inv[5,]
# in regions 193-196:
# within_ancestry_diversity/results/allele_freq/mex2$ zgrep -m 1 "^4" mexicana.allo.withXochi35/*.mafs.gz
regions <- 193:196
#pops <- c("mexicana.allo.withXochi35", paste0("pop", c(18:31,33:35,360:363,365:374)))
pops <- c("mexicana.allo.withXochi35", paste0("pop", c(18:19,21,24:31,34,360:363,365:374))) # add back in #23 later
# add back in ones with problems later: #23 #374 370 366 31 28 21 19
pops <- c("mexicana.allo.withXochi35", paste0("pop", c(18,24:27,29:30,34,360:363,365,367:369,371:373))) 
sites <- do.call(rbind,
                lapply(regions, function(r)
                read.table(paste0("../data/geno_lik/merged_pass1_all_alloMaize4Low_16/allVar_depthFilt/region_", r, ".var.sites"),
                    stringsAsFactors = F, header = F)[,1:2]))
colnames(sites) <- c("chromo", "position")

freqs <- lapply(pops, function(pop) do.call(rbind,
                                            lapply(regions, function(r)
                                              read.table(paste0("results/allele_freq/mex2/", pop, "/region_", r, ".mafs.gz"),
                                                         stringsAsFactors = F, header = T)[ , c(1,2,6,7)])))
names(freqs) <- pops
#colnames(freqs[["mexicana.allo.withXochi35"]]) <- c("chromo", "position", "phat", "nInd")
for (pop in pops){
  colnames(freqs[[pop]]) <- c("chr", "pos", "p", "n")
  freqs[[pop]]$pi <- 2*freqs[[pop]]$p*(1 - freqs[[pop]]$p)
}

# meta data on populations
meta <- unique(read.table("../data/pass1_ids.txt", stringsAsFactors = F, 
                          header = T, sep = "\t")[ , c("popN", "zea", "LOCALITY")])
meta$pop <- paste0("pop",(meta$popN))


pi_within <- data.frame(pop = pops, pi = sapply(pops, function(pop) mean(freqs[[pop]]$pi)))
p_inv4m <- left_join(pi_within, meta, by = "pop") %>%
  ggplot(., aes(x = popN, y = pi, color = LOCALITY)) +
  geom_point()
ggsave(filename = paste0("plots/pi_inv4m.png"), device = "png",
       plot = p_inv4m, 
       width = 10, height = 5, units = "in",
       dpi = 300)
#pi_between <- MAKE PAIRS HERE. CALCULATE FST FOR EACH PAIR OF POPS 



pop30 <- do.call(rbind,
                 lapply(regions, function(r)
                   read.table(paste0("results/allele_freq/mex2/", "pop30", "/region_", r, ".mafs.gz"),
                              stringsAsFactors = F, header = T)[ , c(1,2,6,7)])) %>%
  rename(p.pop30 = phat) %>%
  rename(n.pop30 = nInd)

popMexAllo <- do.call(rbind,
                               lapply(regions, function(r)
                                 read.table(paste0("results/allele_freq/mex2/", "mexicana.allo.withXochi35", "/region_", r, ".mafs.gz"),
                                            stringsAsFactors = F, header = T)[ , c(1,2,6,7)])) %>%
  rename(p.popMexAllo = knownEM) %>% # check this out later
  rename(n.popMexAllo = nInd)
pop367 <- do.call(rbind,
                           lapply(regions, function(r)
                             read.table(paste0("results/allele_freq/mex2/", "pop367", "/region_", r, ".mafs.gz"),
                                        stringsAsFactors = F, header = T)[ , c(1,2,6,7)])) %>%
  rename(p.pop367 = phat) %>%
  rename(n.pop367 = nInd)
all <- dplyr::left_join(sites, pop30, by = c("chromo", "position")) %>%
  dplyr::left_join(., popMexAllo, by = c("chromo", "position")) %>%
  dplyr::left_join(., pop367, by = c("chromo", "position"))
all2 <- do.call(function(d) left_join(., by = c("chromo", "position")),
                list(sites, popMexAllo, pop367))
  
#mexicana.allo.withXochi35
#maize.allo.4Low16

#within_ancestry_diversity/results/allele_freq/maize2$ ln -s ~/hilo/data/geno_lik/merged_pass1_all_alloMaize4Low_16/allVar_depthFilt/maize.allo.4Low16 maize.allo.4Low16
#within_ancestry_diversity/results/allele_freq/mex2$ ln -s ~/hilo/data/geno_lik/merged_pass1_all_alloMaize4Low_16/allVar_depthFilt/mexicana.allo.withXochi35 mexicana.allo.withXochi35

source("rmap_functions.R") # get functions and recomb. map
source("gene_density.R")

dir_in = "../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/input/"
pos_thin <- read.table(paste0(dir_in, "pop18.anc_hmm.input"),
                       stringsAsFactors = F)[ , 1:2] # all pops have same sites
colnames(pos_thin) <- c("chr", "pos")

# calculate map rate for every position in pos_thin
pos_thin$rate <- apply(pos_thin, 1, 
                       function(i) calcMapRate(chr = i["chr"], 
                                               pos = i["pos"]))
write.table(pos_thin, 
            paste0(dir_in, "pos_recomb_rates.txt"), 
            col.names = T, row.names = F, quote = F, sep = "\t")

# calculate gene density for every position in pos_thin
# first by % bp
densCDS = t(apply(pos_thin, 1,
                  function(i) densityCDS(chrom = as.integer(i["chr"]), 
                                         pos = as.integer(i["pos"]))))
write.table(densCDS, 
            paste0(dir_in, "pos_density_CDS.txt"), 
            col.names = T, row.names = F, quote = F, sep = "\t")
  
linkageCDS = t(apply(pos_thin, 1,
                  function(i) linkageCDS(chrom = as.integer(i["chr"]), 
                                         pos = as.integer(i["pos"]))))
write.table(linkageCDS, 
            paste0(dir_in, "pos_linkage_CDS.txt"), 
            col.names = T, row.names = F, quote = F, sep = "\t")


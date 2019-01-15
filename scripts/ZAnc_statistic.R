# this script loads population and individual ancestry data
# and calculates K matrices as well as the ZAnc statistic
# across all sites. Also loads gene density and recomb. rate data
# for all sites with ancestry information.
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(scales)
require(hexbin)
library(reshape2)
library(zoo) # for rolling mean (mean across windows)

# helper file with some useful functions
source("../../covAncestry/forqs_sim/k_matrix.R") # import useful functions

# function for calculating statistic:
make_calcs = function(pop_anc){
  # calculate mean population frequency across all snps
  pop_alpha = apply(pop_anc, 1, mean)
  # calculate K matrix for populations
  pop_K = calcK(ancFreqMatrix = pop_anc, alpha = pop_alpha)
  pop_InvL = calcInvL(pop_K)
  # for each locus, calculate ZAnc statistic
  pop_ZAnc = apply(pop_anc, 2, function(l) ZAnc(ancFreq = l, invL = pop_InvL, alpha = pop_alpha))
  pop_p_ZAnc = calcProbZAnc(ZAnc = pop_ZAnc, nPop = length(pop_alpha), logProb = F)
  return(list(alpha = pop_alpha, K = pop_K, InvL = pop_InvL, ZAnc = pop_ZAnc, p_ZAnc = pop_p_ZAnc))
}

# load needed data:
# get inversions:
inv = read.table("../data/refMaize/inversions/knownInv_v4_coord.txt",
                 stringsAsFactors = F, header = T)
excl_inv_buffer = 1000000 # exclude 1 megabase around inversion
inv$excl_start = inv$start - excl_inv_buffer
inv$excl_end = inv$end + excl_inv_buffer
# which inversions are segregating
# and therefore should be excluded from some analyses?
inv$present = (inv$ID %in% c("Inv4m"))


dir_in = "../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/"
#dir_in = "../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/"
dir_anc = paste0(dir_in, "output_noBoot/anc/")
# color palette for locations (n=13 locations)
location_colors = c("gray10", "deepskyblue1",
  brewer.pal(11, "Spectral"))
# color palette for mex and maize
blues = brewer.pal(n = 9, name = "Blues")[c(6,8)]
yellows = brewer.pal(n = 9, name = "YlOrBr")[c(4,6)]
colors_maize2mex = c(yellows, blues)
labels_maize2mex = c("allopatric_maize", "sympatric_maize", "sympatric_mexicana", "allopatric_mexicana")

# admixture proportions per pop
alphas <- read.table("../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/input/globalAdmixtureByPopN.txt",
                     stringsAsFactors = F, header = F)
colnames(alphas) <- c("popN", "alpha_maize", "alpha_mex")
# metadata for each individual, including pop association
pop_elev <- read.table("../data/riplasm/gps_and_elevation_for_sample_sites.txt",
                       stringsAsFactors = F, header = T, sep = "\t") %>%
  dplyr::select(., popN, ELEVATION)
meta <- read.table("../data/pass1_ids.txt", stringsAsFactors = F, 
                   header = T, sep = "\t") %>%
  filter(., est_coverage >= 0.05) %>%
  left_join(., alphas, by = c("popN")) %>%
  left_join(., pop_elev, by = c("popN"))

# individual by individual K matrix
included_inds = meta %>%
  filter(alpha_maize > 0 & alpha_mex > 0) %>%
  .[order(.$popN), ]
pops = unique(included_inds[order(included_inds$popN), c("popN", "zea", "LOCALITY", "alpha_maize", "alpha_mex", "ELEVATION")])
pops$zea_SHORT = abbreviate(pops$zea, min = 4, use.classes = F)
pops$LOC_SHORT = abbreviate(pops$LOCALITY, min = 6, method = "both.sides", use.classes = F)
#write.table(pops$popN, paste0(dir_in, "/included_pops.txt"),
#            col.names = F, row.names = F, quote = F)

# read in population ancestry input files from calc_genomewide_pop_anc_freq.R
pop_anc_list = lapply(pops$popN, function(pop) read.table(paste0(dir_anc, "/pop", pop, ".anc.freq"), 
                                                stringsAsFactors = F))

# combine population ancestry frequencies into a matrix
# where rows are populations and columns are snps
all_anc = t(do.call(cbind, pop_anc_list))
rownames(all_anc) <- paste(pops$zea_SHORT, pops$LOC_SHORT, sep = ".")

# run K-matrix and ZAnc calculations:
all = make_calcs(pop_anc = all_anc)
# just maize separately:
maize_anc = all_anc[pops$zea == "maize", ]
maize = make_calcs(maize_anc)
# just mexicana separately:
mex_anc = all_anc[pops$zea == "mexicana", ]
mex = make_calcs(mex_anc)


# read in individual ancestry input files from calc_genomewide_pop_anc_freq.R
ind_anc_list = lapply(pops$popN, function(pop) read.table(paste0(dir_anc, "/pop", pop, ".anc.ind"), 
                                                          stringsAsFactors = F))
# How many haplotypes were sequenced per population?
pops$N_haps = unlist(lapply(ind_anc_list, 
                            function(i) dim(i)[2]*2))

# combine population ancestry frequencies into a matrix
# where rows are populations and columns are snps
all_ind_anc = t(do.call(cbind, ind_anc_list))
#rownames(all_ind_anc) <- paste(substr(included_inds$zea, 1, 4), substr(included_inds$LOCALITY, 1, 2), sep = ".")
rownames(all_ind_anc) <- included_inds$ID

# calculate individual K matrix
all_ind = make_calcs(pop_anc = all_ind_anc)

# script rmap.R generates recombination rates for each position in thinned positions
#pos <- read.table(
#            paste0(dir_in, "input/pos_recomb_rates.txt"), 
#            header = T, stringsAsFactors = F)
# use if no recombination rates have been calculated
#pos = read.table(paste0(dir_in, "/output_noBoot/HILO1.posterior"),
#                 stringsAsFactors = F, header = T)[ , 1:2] %>%
#  rename(pos = position)

# I'll get position info from var.sites used for ancestry_hmm
# and then gene density from the 0.1cM windows around each SNP
sites <- do.call(rbind, 
               lapply(1:10, function(i)
                 read.table(paste0(dir_in, "../chr", i, ".var.sites"),
                            header = F, stringsAsFactors = F)))
colnames(sites) <- c("chr", "pos", "major", "minor")
wind <- read.table(paste0(dir_in, "/../windows0.1cM/gene_overlap.bed"),
                 header = F, stringsAsFactors = F)
colnames(wind) <- c("chr_wind", "start", "end", 
                    "width_cM", "coding_bp",
                    "n_CDS", "width_bp", # note number CDS is number contigous coding sequences (may be more than one gene (?))
                    "perc_coding")
pos <- cbind(sites, wind)
# calculate local recombination rate
# in cM/Mb
pos$rate <- pos$width_cM/(pos$width_bp/10^6)
hist(pos$rate)
summary(pos$rate)
# quantiles for defining 'low' and 'high' recomb.
low_quant = .2
high_quant = .8
quantile(pos$rate, c(low_quant, high_quant)) # 0.2 and 4.0

# low recombination rate regions <= low quantile:
low_r_anc = all_anc[ , pos$rate <= quantile(pos$rate, low_quant)]
low_r = make_calcs(low_r_anc)
low_r_ind_anc = all_ind_anc[ , pos$rate <= quantile(pos$rate, low_quant)]
low_r_ind = make_calcs(low_r_ind_anc)
# high recombination rate regions >= high quantile:
high_r_anc = all_anc[ , pos$rate >= quantile(pos$rate, high_quant)]
high_r = make_calcs(high_r_anc)
high_r_ind_anc = all_ind_anc[ , pos$rate >= quantile(pos$rate, high_quant)]
high_r_ind = make_calcs(high_r_ind_anc)

# which positions / loci are within present inversions?
# for each inversion present, create a column in pos table
# for whether that SNP is within the inversion
for (i in unique(inv[inv$present, "ID"])){
  pos[ , paste0("inv_", i)] <- (pos$chr == inv[inv$ID == i, "chrom"] &
                                  pos$pos >= inv[inv$ID == i, "excl_start"] &
                                  pos$pos <= inv[inv$ID == i, "excl_end"])
}
# is this position within any segregating inversions?
pos$inv_any <- dplyr::select(pos, starts_with("inv_")) %>%
  apply(., 1, any) # across all the inversions, are any values TRUE?

# put all the K and ZAnc calculations together in one list

# individual and population ancestry calls
ancestry = list(ind = all_ind_anc, pop = all_anc)

# first with all loci, including those within inversions
with_inv = list(ind = list(high_r = high_r_ind, low_r = low_r_ind, all = all_ind), 
                pop = list(high_r = high_r, low_r = low_r, all = all))

# then do the same for positions not within any known inversions
no_inv = lapply(ancestry, function(i)
  list(high_r = i[ , !pos$inv_any & pos$rate >= quantile(pos$rate, high_quant)], 
       low_r = i[ , !pos$inv_any & pos$rate <= quantile(pos$rate, low_quant)], 
       all = i[ , !pos$inv_any ]) %>% 
    lapply(., make_calcs)) 


#note:
# I am not sure whether this is the appropriate
# alpha to recalculate it for just high or low recombining regions --
# I think I should be using a global alpha from genomewide
# or maybe weighing it differently than just taking mean across thinned snps
# because SNPs kept are more likely in low recomb. areas
# BUT we do think the alphas from low recombination are
# being affected by selection; in finding e.g. mexicana
# outlier loci in maize it's conservative to take the genomewide
# avg. in regions of low recombination but maybe anti-conservative 
# for potenital outliers in regions of high recombination




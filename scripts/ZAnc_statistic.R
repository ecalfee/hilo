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

# get inversions:
inv = read.table("../data/refMaize/inversions/knownInv_v4_coord.txt",
                 stringsAsFactors = F, header = T)
excl_inv_buffer = 1000000 # exclude 1 megabase around inversion

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
all = make_calcs(pop_anc = all_anc)
# just maize separately:
maize_anc = all_anc[pops$zea == "maize", ]
maize = make_calcs(maize_anc)
# just mexicana separately:
mex_anc = all_anc[pops$zea == "mexicana", ]
mex = make_calcs(mex_anc)

# plot population K as a correlation matrix
png(paste("../plots/K_corrplot_combined.png"), # saves plot as pin in ../plots/
    height = 7, width = 7, units = "in", res = 150)
corrplot(cov2cor(all$K), method="shade", main = "mex/maize pop anc corr.")
dev.off()

# plot population K raw covariance matrix
png(paste("../plots/K_covplot_combined.png"), # saves plot as pin in ../plots/
    height = 7, width = 7, units = "in", res = 150)
#heatmap(all$K, Colv = NA, Rowv = NA, main = "mex/maize pop anc cov.",
#        col= colorRampPalette(brewer.pal(8, "Blues"))(25))
col_pop_location=location_colors[my_group=as.numeric(as.factor(pops$LOCALITY))]
heatmap(all$K, Colv = NA, Rowv = NA, main = "mex/maize pop anc cov.",
        col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
        RowSideColors = col_pop_location,
        scale = "none")
dev.off()

# read in individual ancestry input files from calc_genomewide_pop_anc_freq.R
ind_anc_list = lapply(pops$popN, function(pop) read.table(paste0(dir_anc, "/pop", pop, ".anc.ind"), 
                                                          stringsAsFactors = F))

# combine population ancestry frequencies into a matrix
# where rows are populations and columns are snps
all_ind_anc = t(do.call(cbind, ind_anc_list))
#rownames(all_ind_anc) <- paste(substr(included_inds$zea, 1, 4), substr(included_inds$LOCALITY, 1, 2), sep = ".")
rownames(all_ind_anc) <- included_inds$ID

# calculate individual K matrix
all_ind = make_calcs(pop_anc = all_ind_anc)

# plot individual K, where individuals are grouped (for visualization) by pop,
# as a correlation matrix
png(paste("../plots/K_corrplot_ind_combined.png"), # saves plot as pin in ../plots/
    height = 21, width = 21, units = "in", res = 150)
corrplot(cov2cor(all_ind$K), method="shade", main = "mex/maize ind anc corr.")
dev.off()

# plot individual K as a raw covariance matrix
png(paste("../plots/K_covplot_ind_combined.png"), # saves plot as pin in ../plots/
    height = 21, width = 21, units = "in", res = 150)
col_ind_location=location_colors[my_group = 
                                                 as.numeric(as.factor(included_inds$LOCALITY))]
heatmap(all_ind$K, 
        Colv = NA, Rowv = NA, 
        main = "mex/maize ind anc cov.",
        col = colorRampPalette(brewer.pal(8, "Blues"))(25),
        RowSideColors = col_ind_location,
        scale = "none")
dev.off()



# I don't have a good idea yet how much of the ancestry covariance across all maize has to do with selection
# at the same loci for mexicana introgression.
# comparing K matrices for low and high recombining regions of the genome
# will help answer this question of selection vs. demography

# A first step is to get recombination rates for every thinned position
# and plot mean ancestry by recombination rate 
# and visualize K by high/low recomb. where alpha is a genomewide mean 
# and visualize K by high/low recomb. where alpha is a recombination bin mean
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

# low vs. high recomb. population level K 
png(paste("../plots/K_covplot_pop_combined_low_recomb_rate.png"), # saves plot as pin in ../plots/
    height = 7, width = 7, units = "in", res = 150)
heatmap(low_r$K, Colv = NA, Rowv = NA, main = "low recomb - pop anc cov.",
        col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
        RowSideColors = col_pop_location,
        scale = "none")
dev.off()

png(paste("../plots/K_covplot_pop_combined_high_recomb_rate.png"), # saves plot as pin in ../plots/
    height = 7, width = 7, units = "in", res = 150)
heatmap(high_r$K, Colv = NA, Rowv = NA, main = "high recomb - pop anc cov.",
        col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
        RowSideColors = col_pop_location,
        scale = "none")
dev.off()


# low vs. high recomb. ind. level K
png(paste("../plots/K_covplot_ind_combined_low_recomb_rate.png"), # saves plot as pin in ../plots/
    height = 14, width = 14, units = "in", res = 150)
heatmap(low_r_ind$K, Colv = NA, Rowv = NA, main = "low recomb - ind anc cov.",
        col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
        RowSideColors = col_ind_location,
        scale = "none")
dev.off()

png(paste("../plots/K_covplot_ind_combined_high_recomb_rate.png"), # saves plot as pin in ../plots/
    height = 14, width = 14, units = "in", res = 150)
heatmap(high_r_ind$K, Colv = NA, Rowv = NA, main = "high recomb - ind anc cov.",
        col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
        RowSideColors = col_ind_location,
        scale = "none")
dev.off()

# correlation matrices low vs. high recomb.
# population level
png(paste("../plots/K_corrplot_pop_combined_low_recomb_rate.png"),
    height = 7, width = 7, units = "in", res = 150)
corrplot(cov2cor(low_r$K), method="shade", main = "Low recomb. pop anc corr.")
dev.off()
png(paste("../plots/K_corrplot_pop_combined_high_recomb_rate.png"),
    height = 7, width = 7, units = "in", res = 150)
corrplot(cov2cor(high_r$K), method="shade", main = "High recomb. pop anc corr.")
dev.off()



# individual level
png(paste("../plots/K_corrplot_ind_combined_low_recomb_rate.png"),
    height = 14, width = 14, units = "in", res = 150)
corrplot(cov2cor(low_r_ind$K), method="shade", main = "Low recomb. ind anc corr.")
dev.off()
png(paste("../plots/K_corrplot_ind_combined_high_recomb_rate.png"),
    height = 14, width = 14, units = "in", res = 150)
corrplot(cov2cor(high_r_ind$K), method="shade", main = "High recomb. ind anc corr.")
dev.off()

# I can plot the covariances from the K matrices 
# as box plots to compare covariances
# for specific pairs of pops

# do pairs of maize-mex populations at the same location 
# covary more with each other than random pairs?
# which populations have mex-maize pairs?
pop_pairs = pops$LOCALITY[duplicated(pops$LOCALITY)]
pops$paired = pops$LOCALITY %in% pop_pairs
pop_pairs_SHORT = unique(pops$LOC_SHORT[pops$LOCALITY %in% pop_pairs])
# How many haplotypes were sequenced per population
# (so I can do a finite sample size correction for
# the variance within pop to fairly compare it with
# covariances between populations)
pops$N_haps = unlist(lapply(ind_anc_list, 
                            function(i) dim(i)[2]))


# flatten or 'melt' covariance matrices K
kmelt <- function(K, haps_per_pop, regions = NULL){ 
  # haps per pop is to do the finite sample size correction
  # for the variance (on the diagonal) so it can be compared
  # to covariances between populations
  
  # not quite : want matrix 1's everywhere else but N/N-1 on diagonals--
  K_corrected = K*(matrix(1, dim(K)[1], dim(K)[2]) + diag(haps_per_pop))
  k_melted = reshape2::melt(K_corrected, 
                         value.name = "covariance")[upper.tri(K_corrected, diag = T),] %>%
    separate(., Var1, c("zea_SHORT_1", "LOC_SHORT_1")) %>%
    separate(., Var2, c("zea_SHORT_2", "LOC_SHORT_2")) %>%
    mutate(., same_location = LOC_SHORT_1 == LOC_SHORT_2) %>%
    mutate(., mex_vs_maize = zea_SHORT_1 != zea_SHORT_2) %>%
    mutate(., have_pairs = LOC_SHORT_1 %in% pop_pairs_SHORT & LOC_SHORT_2 %in% pop_pairs_SHORT)
  if (! is.null(regions)){ # label these results
    k_melted$regions <- regions
  }
  return(k_melted)
  }

melted_ks <- do.call(rbind,
             lapply(c("all", "low_r", "high_r"),
                    function(i) kmelt(get(i)$K, regions = i)))
# show patterns for covarience between maize-mex pairs
melted_ks %>%
  filter(., have_pairs & mex_vs_maize) %>%
  ggplot(., aes(y = covariance, x = same_location,
                      fill = regions)) + 
  geom_boxplot() +
    ggtitle("covariance mex-maize pairs at same and diff locations")
# show patterns of covariance for maize-maize pairs
melted_ks %>%
  filter(., (!mex_vs_maize) & zea_SHORT_1 == "maiz") %>%
  ggplot(., aes(y = covariance, x = same_location,
                      fill = regions)) + 
  geom_boxplot() +
    ggtitle("covariance maize-maize pairs at same and diff locations")
# show patterns of covariance for mex-mex pairs
melted_ks %>%
  filter(., (!mex_vs_maize) & zea_SHORT_1 == "mexi") %>%
  ggplot(., aes(y = covariance, x = same_location,
                      fill = regions)) + 
  geom_boxplot() +
    ggtitle("covariance mex-mex pairs at same and diff locations")


# pairs_k is a k matrix with covariances between mexicana pops (rows)
# and maize pops (cols), where diagonal is same location,
# and only pops with a pair are included in the matrix
for (i in c("all", "low_r", "high_r")){
    k = get(i)
  pairs_k_unordered = k$K[pops$zea == "mexicana" & pops$paired, 
                            pops$zea == "maize" & pops$paired]
  k$pairs_k = pairs_k_unordered[order(rownames(pairs_k_unordered)), 
                        order(colnames(pairs_k_unordered))]
  # make dataframe with columns 'paired', 'unpaired', 'nonpair'
  cov_pairs = diag(k$pairs_k)
  cov_nonpairs = pairs_k[col(k$pairs_k) != row(k$pairs_k)]
  boxplot(cov_pairs, cov_nonpairs,
          main = i)
}






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




# plot estimated mexicana ancestry
# by high vs. low recomb. rate regions
png(paste("../plots/alpha_by_recomb_rate.png"), # saves plot as pin in ../plots/
    height = 5, width = 7, units = "in", res = 150)
plot(y = low_r$alpha, x = high_r$alpha,
     main = "mean pop. mexicana ancestry by recomb. rate",
     ylab = "alpha low r",
     xlab = "alpha high r",
     col = ifelse(pops$zea == "mexicana", "blue", "orange"))
legend(x = "topleft",
  legend = c("maize", "mexicana", "1-to-1 line"), 
       col = c("orange", "blue", "black"),
       pch = c(1, 1, 5))
abline(a = 0, b = 1)
dev.off()





# how much of the genome appears to be under selection?
sum(maize[["p_ZAnc"]] < .01)/sum(maize[["ZAnc"]] > 0)
# ~ 9% of teosinte-biased ancestry segments appear to be selected for mexicana upon rough approximation
maize_bound = qnorm(p = .01, mean = 0, sd = sqrt(nrow(maize_anc)), lower.tail = F)
sum(abs(maize[["ZAnc"]]) > maize_bound)/ncol(maize_anc) # 3.6% genome overall appears to be under selection
sum(maize[["ZAnc"]] > maize_bound)/ncol(maize_anc) # nearly all of it in the teosinte direction
mex_bound = qnorm(p = .01, mean = 0, sd = sqrt(nrow(mex_anc)), lower.tail = F)
all_bound = qnorm(p = .01, mean = 0, sd = sqrt(nrow(all_anc)), lower.tail = F)
sum(mex[["ZAnc"]] > mex_bound)/ncol(mex_anc) # hitting boundary?
sum(mex[["ZAnc"]] < -mex_bound)/ncol(mex_anc) # 1-2% maize biased


# histogram of ZAnc statistic
pop_results = list(all, maize, mex)
pop_name = list("all", "maize", "mex")
bounds = list(all_bound, maize_bound, mex_bound)
for (i in 1:length(pop_results)){
  pr <- pop_results[[i]]
  name <- pop_name[[i]]
  png(paste("../plots/ZAnc_hist_", name , ".png"), # saves plot as pin in ../plots/
      height = 5, width = 5, units = "in", res = 150)
  hist(pr[["ZAnc"]], main = paste("ZAnc", name))
  mtext(paste0("p<.01 outliers maize:", round(sum(pr[["ZAnc"]] < -bounds[[i]])/ncol(all_anc), 3), 
               " mex:", round(sum(pr[["ZAnc"]] > bounds[[i]])/ncol(all_anc), 3))) # I suspect outliers are driven mainly by maize
  abline(v = c(bounds[[i]], -bounds[[i]]), col = "red")
  dev.off()
}
# plot p-values for maize ZAnc
plot(maize[["ZAnc"]], maize[["p_ZAnc"]], cex = .1, 
     col = ifelse(maize[["p_ZAnc"]] < .01 | maize[["p_ZAnc"]] > .99, "red", "blue"),
     main = "p_values for ZAnc statistic across sites in maize")


# plot alpha NGSAdmix vs. alpha from ancestry_hmm
pops$alpha_ancestry_hmm = all[["alpha"]]
ggplot(pops, aes(x = alpha_mex, y = alpha_ancestry_hmm)) +
  geom_point(aes(col = LOCALITY)) +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  ggtitle("genome-wide mexicana ancestry: NGSadmix vs. ancestry_hmm")
ggsave("../plots/alpha_est_NGSadmix_vs_ancestry_hmm.png", height = 5, width = 5, units = "in", device = "png")

# now using positions I plot ZAnc by chromosome:

d = data.frame(pos, ZAnc_maize = maize[["ZAnc"]], ZAnc_mex = mex[["ZAnc"]], ZAnc_all = all[["ZAnc"]])
d %>%
  #filter(., chr == "4") %>%
  ggplot(., aes(x = pos, y = ZAnc_maize)) + 
  geom_point(cex= .1, alpha = .5, col = "orange") + 
  facet_wrap(~chr)
ggsave("../plots/ZAnc_maize_whole_genome.png", height = 15, width = 15, units = "in", device = "png")
d %>%
  #filter(., chr == "4") %>%
  ggplot(., aes(x = pos, y = ZAnc_mex)) + 
  geom_point(cex= .1, alpha = .5, col = "blue") + facet_wrap(~chr)
ggsave("../plots/ZAnc_mex_whole_genome.png", height = 15, width = 15, units = "in", device = "png")
d %>%
  #filter(., chr == "4") %>%
  ggplot(., aes(x = pos, y = ZAnc_all)) + 
  geom_point(cex= .1, alpha = .5, col = "black") + facet_wrap(~chr)
ggsave("../plots/ZAnc_combined_whole_genome.png", height = 15, width = 15, units = "in", device = "png")


# plot ZAnc against recombination rates in 0.1cM windows around SNPs
d %>%
  filter(., chr == "4") %>%
  ggplot(., aes(x = pos, y = rate)) + 
  geom_point(cex= .1, alpha = .5, col = "black") +
  geom_smooth(col = "blue")
d %>%
  #filter(., chr == "4") %>%
  ggplot(., aes(x = rate, y = ZAnc_maize)) + 
  geom_point(cex= .1, alpha = .5, col = "black") +
  geom_smooth(col = "blue")

# group rates by quantiles for barplots:
d$bin_rate = cut(d$rate, # include full range of depth in cut
                 breaks = quantile(d$rate, p = seq(0, 1, by = .2)),
                 right = T,
                 include.lowest = T)
d$coding_density = d$coding_bp/d$width_cM
d$bin_coding_density = cut(d$coding_density,
                           breaks = quantile(d$coding_density, p = seq(0, 1, by = .2)),
                           right = T,
                           include.lowest = T)
d$meanAlpha_maize = colMeans(maize_anc)
d$meanAlpha_mex = colMeans(mex_anc)
# exclude inversion area
d$excl_inv <- d$chr == inv[inv$ID=="Inv4m", "chrom"] & 
  d$pos > inv[inv$ID == "Inv4m", "start"] - excl_inv_buffer &
  d$pos < inv[inv$ID == "Inv4m", "end"] + excl_inv_buffer
d_split <- d %>%
  filter(., !excl_inv) %>%
  gather(., "group_alpha", "alpha", 
         c("meanAlpha_mex", "meanAlpha_maize")) %>%
  gather(., "group_ZAnc", "ZAnc",
         c("ZAnc_mex", "ZAnc_maize", "ZAnc_all"))

# make bar plot
box_rate_ZAnc <- 
  ggplot(d_split, aes(y = ZAnc, x = bin_rate,
                fill = group_ZAnc)) + 
  scale_fill_manual(values = c("grey", colors_maize2mex[c(1,4)])) +
  geom_boxplot() +
  facet_wrap(~group_ZAnc) +
  labs(x = "recombination rate (0.1cM winds)", 
       y = "ZAnc statistic by group",
       main = "Effect of local recomb. on ZAnc statistic")
plot(box_rate_ZAnc)
ggsave("../plots/box_rate_ZAnc_0.1cM_windows.png", 
       height = 15, width = 20, 
       units = "in", device = "png")

# make bar plot for mean ancestry mexicna
box_rate_alpha <- 
  ggplot(d_split, aes(y = alpha, x = bin_rate,
                      fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_boxplot() +
  facet_wrap(~group_alpha) +
  labs(x = "recombination rate (0.1cM winds)", 
       y = "Mean mex. ancestry by group",
       main = "Effect of local recomb. on mean mex ancestry (alpha)")
plot(box_rate_alpha)
ggsave("../plots/box_rate_alpha_0.1cM_windows.png", 
       height = 15, width = 20, 
       units = "in", device = "png")

# now for gene density: make bar plot for mean ancestry mexicna
box_coding_density_alpha <- 
  ggplot(d_split, aes(y = alpha, x = bin_coding_density,
                      fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_boxplot() +
  facet_wrap(~group_alpha) +
  labs(x = "coding density (bp/cM)", 
       y = "Mean mex. ancestry by group",
       main = "Effect of local coding desnity on mean mex ancestry (alpha)")
plot(box_coding_density_alpha)
ggsave("../plots/box_coding_density_alpha_0.1cM_windows.png", 
       height = 15, width = 20, 
       units = "in", device = "png")


# plot ZAnc against gene density in 0.1cM windows around SNPs
d %>%
  #filter(., chr == "4") %>%
  ggplot(., aes(x = perc_coding, y = ZAnc_maize)) + 
  geom_point(cex= .1, alpha = .5, col = "black") +
  geom_smooth(col = "blue")


# plot recomb. rate against gene density in 0.1cM windows around SNPs
d %>%
  #filter(., chr == "4") %>%
  ggplot(., aes(x = rate, y = perc_coding)) + 
  geom_point(cex= .1, alpha = .5, col = "black") +
  geom_smooth(col = "blue")

# plot ancestry by recombination rate
#rownames(d)<-NULL
# add in ancestry information for each population 
# at each site
d_anc <- bind_cols(d, data.frame(t(maize_anc)), 
                   data.frame(t(mex_anc))) %>%
  gather(., "pop", "alpha_mex", 17:37)
# for individual pops, plot recomb. rate vs. ancestry
d_anc %>%
  filter(., chr == "4") %>%
  filter(., pop == "maiz.Na" | pop == "maiz.Op") %>%
  ggplot(., aes(x = rate, y = alpha_mex)) +
  geom_point(cex= .1, alpha = .5, col = "black") +
  geom_smooth(col = "blue") +
  facet_wrap(~pop)
ggsave("../plots/anc_by_recomb_rate_0.1cM_windows.png", 
       height = 20, width = 25, 
       units = "in", device = "png")
# and plot coding density vs. ancestry
d_anc %>%
  #filter(., chr == "4") %>%
  #filter(., pop == "maiz.Na" | pop == "maiz.Op") %>%
  ggplot(., aes(x = coding_bp, y = alpha_mex)) +
  geom_point(cex= .1, alpha = .5, 
             col = "black") +
  geom_smooth(col = "blue") +
  facet_wrap(~pop)
ggsave("../plots/anc_by_coding_bp.png", 
       height = 20, width = 25, 
       units = "in", device = "png")
# percent coding bp
d_anc %>%
  #filter(., chr == "4") %>%
  #filter(., pop == "maiz.Na" | pop == "maiz.Op") %>%
  ggplot(., aes(x = coding_bp/width_cM, y = alpha_mex)) +
  geom_bin2d() +
  #geom_point(cex= .1, alpha = .5, 
  #           col = "black") +
  geom_smooth(col = "blue") +
  facet_wrap(~pop)
ggsave("../plots/anc_by_coding_bp_per_cM.png", 
       height = 20, width = 25, 
       units = "in", device = "png")

# plain linked bp within 0.1cM windows
p_bp <- d_anc %>%
  #filter(., chr == "4") %>%
  #filter(., pop == "maiz.Na" | pop == "maiz.Op") %>%
  ggplot(., aes(x = width_bp, y = alpha_mex)) 
#p_bp + geom_hex() # try hexbin out with smaller data
p_bp + geom_hex() +
  geom_smooth(col = "blue") +
  facet_wrap(~pop)
p_bp + geom_bin2d() +
  geom_smooth(col = "blue") +
  facet_wrap(~pop)
#p_bp + geom_point(cex= .1, alpha = .05, 
#             col = "black") +
#  geom_smooth(col = "blue") +
#  facet_wrap(~pop)
ggsave("../plots/anc_by_width_bp_0.1cM_windows_bins.png", 
       height = 20, width = 25, 
       units = "in", device = "png")

# Now get smoothed ZAnc statistic plots 
# (taking mean across windows of 50 thinned positons ~ 50-100 kb)
# geom_smooth (method = "gam" (loess is too memory intense for > 1000 points))
# geom_smooth uses a regression method (of your choice) ..
# to actually just use raw mean of a window I'll use roll_mean() from package zoo
w = 100 # window size (# thinned snps)
d %>%
  #filter(., chr == "4") %>%
  # calculate rolling (or windowed) mean ZAnc
  mutate(., ZAnc_maize_wind = zoo::rollmean(ZAnc_maize, w, fill = list(NA, NULL, NA))) %>%
  .[c(T, rep(F, as.integer(w/5) - 1)), ] %>% # plot every 5th point
  ggplot(., aes(x = pos, y = ZAnc_maize)) + 
  geom_line(size = .5, cex= .1, alpha = .5, col = "orange") + 
  #geom_line(size = 1.2, cex= .1, alpha = .5, col = "orange") +
  facet_wrap(~chr)

# compare individual to population results for 
# ZAnc
# first in maize
# just maize separately:
maize_ind_anc = all_ind_anc[included_inds$zea == "maize", ]
maize_ind = make_calcs(maize_ind_anc)
# just mexicana separately:
mex_ind_anc = all_ind_anc[included_inds$zea == "mexicana", ]
mex_ind = make_calcs(mex_ind_anc)
# add individual ZAnc estimates to d
d1 <- d %>%
  mutate(., ZAnc_maize_ind = maize_ind$ZAnc) %>%
  mutate(., ZAnc_mex_ind = mex_ind$ZAnc) %>%
  mutate(., ZAnc_all_ind = all_ind$ZAnc)%>%
  mutate(., maize_pop_anc_freq = apply(maize_anc, 2, mean)) %>%
  mutate(., mex_pop_anc_freq = apply(mex_anc, 2, mean)) %>%
  mutate(., all_pop_anc_freq = apply(all_anc, 2, mean)) %>%
  mutate(., maize_ind_anc_freq = apply(maize_ind_anc, 2, mean)) %>%
  mutate(., mex_ind_anc_freq = apply(mex_ind_anc, 2, mean)) %>%
  mutate(., all_ind_anc_freq = apply(all_ind_anc, 2, mean))
# plot ZAnc for maize from individual-by-individual K
d1 %>%
  #filter(., chr == "4") %>%
  # calculate rolling (or windowed) mean ZAnc
  mutate(., ZAnc_maize_ind_wind = zoo::rollmean(ZAnc_maize_ind, w, fill = list(NA, NULL, NA))) %>%
  .[c(T, rep(F, as.integer(w/5) - 1)), ] %>% # plot every 5th point
  ggplot(., aes(x = pos, y = ZAnc_maize_ind)) + 
  geom_line(size = .5, alpha = .5, col = "coral") + 
  facet_wrap(~chr)
ggsave("../plots/ZAnc_maize_whole_genome.png", height = 15, width = 15, units = "in", device = "png")
# mex individual-by-individual ZAnc
d1 %>%
  #filter(., chr == "4") %>%
  # calculate rolling (or windowed) mean ZAnc
  mutate(., ZAnc_mex_ind_wind = zoo::rollmean(ZAnc_mex_ind, w, fill = list(NA, NULL, NA))) %>%
  .[c(T, rep(F, as.integer(w/5) - 1)), ] %>% # plot every 5th point
  ggplot(., aes(x = pos, y = ZAnc_mex_ind)) + 
  geom_line(size = .5, alpha = .5, col = "cornflowerblue") + 
  facet_wrap(~chr)
ggsave("../plots/ZAnc_mex_ind_whole_genome.png", height = 15, width = 15, units = "in", device = "png")

# plot raw ancestry frequencies (across individuals)
d1 %>%
  #filter(., chr == "4") %>%
  # calculate rolling (or windowed) mean anc freq
  mutate(., maize_ind_anc_freq_wind = zoo::rollmean(maize_ind_anc_freq, w, fill = list(NA, NULL, NA))) %>%
  .[c(T, rep(F, as.integer(w/5) - 1)), ] %>% # plot every 5th point
  ggplot(., aes(x = pos, y = maize_ind_anc_freq_wind)) + 
  geom_line(size = .5, alpha = .5, col = "red") + 
  facet_wrap(~chr)
ggsave("../plots/ZAnc_maize_ind_anc_freq_whole_genome.png", height = 15, width = 15, units = "in", device = "png")
# averaging population means for ancestry
d1 %>%
  #filter(., chr == "4") %>%
  # calculate rolling (or windowed) mean anc freq
  mutate(., maize_pop_anc_freq_wind = zoo::rollmean(maize_pop_anc_freq, w, fill = list(NA, NULL, NA))) %>%
  .[c(T, rep(F, as.integer(w/5) - 1)), ] %>% # plot every 5th point
  ggplot(., aes(x = pos, y = maize_pop_anc_freq_wind)) + 
  geom_line(size = .5, alpha = .5, col = "darkred") + 
  facet_wrap(~chr)
ggsave("../plots/ZAnc_maize_pop_anc_freq_whole_genome.png", height = 15, width = 15, units = "in", device = "png")

# I need to plot some ZAnc peaks & raw inferred ancestry tracts

# mex individual-by-individual ZAnc
d1 %>%
  #filter(., chr == "4") %>%
  # calculate rolling (or windowed) mean ZAnc
  mutate(., ZAnc_mex_ind_wind = zoo::rollmean(ZAnc_mex_ind, w, fill = list(NA, NULL, NA))) %>%
  .[c(T, rep(F, as.integer(w/5) - 1)), ] %>% # plot every 5th point
  ggplot(., aes(x = pos, y = ZAnc_mex_ind)) + 
  geom_line(size = .5, alpha = .5, col = "navy") + 
  facet_wrap(~chr)
ggsave("../plots/ZAnc_mex_ind_whole_genome.png", height = 15, width = 15, units = "in", device = "png")




# plot individual and population ZAnc together for comparison
# I should also add plotting just raw mean ancestry freq across pops
input366 <- read.table(paste0(dir_in, "input/pop366.anc_hmm.input"),
                       stringsAsFactors = F)
# large popuation (9 individuals)
colnames(input366) <- c("chr", "pos",
                        "maize_A", "maize_a",
                        "mex_A", "mex_a",
                        "rmap_pos", 
                        paste(paste("ind", sapply(1:9, function(x) rep(x, 2)), sep = "_"), c("A", "a"), sep = "_"))
input366 <- input366 %>% # get coverage counts for ref and tot pop
  mutate(., maize_tot = maize_A + maize_a) %>%
  mutate(., mex_tot = mex_A + mex_a) %>%
  mutate(., admix_tot = apply(.[ , 8:25], 1, sum))
# Sanity checks:
h1 <- input366 %>%
  gather(., "ref", "depth", c("mex_tot", "maize_tot")) %>%
  ggplot(., aes(x = depth, fill = ref)) +
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_histogram(position = "identity", alpha = .75) #+ facet_wrap(~ref)
ggsave("../plots/hist_mex_maize_ref_coverage_thinned_pos.png", 
       h1,
       height = 7, width = 10, 
       units = "in", device = "png")


png(paste("../plots/hist_maize_ref_coverage.png"),
    height = 5, width = 5, units = "in", res = 150)
hist(input366$maize_tot, main = "Hist total coverage maize ref")
dev.off()
png(paste("../plots/hist_mex_ref_coverage.png"),
    height = 5, width = 5, units = "in", res = 150)
hist(input366$mex_tot, main = "Hist total coverage mex ref")
dev.off()
png(paste("../plots/ZAnc_by_mex_ref_coverage.png"),
    height = 5, width = 5, units = "in", res = 150)
plot(x = maize$ZAnc, y = input366$mex_tot,
     main = "ZAnc by mex ref. depth",
     col = scales::alpha("orange", .1))
dev.off()
png(paste("../plots/ZAnc_by_maize_ref_coverage.png"),
    height = 5, width = 5, units = "in", res = 150)
plot(x = maize$ZAnc, y = input366$maize_tot, 
     main = "ZAnc by maize ref. depth",
     col = scales::alpha("orange", .1))
dev.off()
png(paste("../plots/ZAnc_by_pop366_tot_coverage_OUTLIERS.png"),
    height = 5, width = 5, units = "in", res = 150)
plot(x = maize$ZAnc, y = input366$admix_tot, 
     main = "ZAnc by total depth pop 366",
     col = scales::alpha("orange", .1))
dev.off()
png(paste("../plots/ZAnc_by_pop366_tot_coverage_NO_OUTLIERS.png"),
    height = 5, width = 5, units = "in", res = 150)
plot(x = maize$ZAnc, y = input366$admix_tot, 
     main = "ZAnc by total depth pop 366",
     ylim = c(0,50),
     col = scales::alpha("orange", .1))
dev.off()

png(paste("../plots/ZAnc_mex_by_maize_ref_coverage.png"),
    height = 5, width = 5, units = "in", res = 150)
plot(x = mex$ZAnc, y = input366$maize_tot, 
     main = "ZAnc by maize ref. depth",
     col = scales::alpha("blue", .1))
dev.off()
png(paste("../plots/ZAnc_mex_by_mex_ref_coverage.png"),
    height = 5, width = 5, units = "in", res = 150)
plot(x = mex$ZAnc, y = input366$mex_tot, 
     main = "ZAnc by mex ref. depth",
     col = scales::alpha("blue", .1))
dev.off()
# better to bin by depth and then layer on ZAnc
bin_breaks <- seq(0, 40, by = 2)
ref_depth <- data.frame(chr = pos$chr,
                        pos = pos$pos,
                        mex_ref = input366$mex_tot,
                        maize_ref = input366$maize_tot,
                        mex = mex$ZAnc,
                        maize = maize$ZAnc) %>%
  gather(., "ref", "tot_depth", c("mex_ref", "maize_ref")) %>%
  mutate(tot_depth_bin = cut(tot_depth, # include full range of depth in cut
                   breaks = bin_breaks,
                   right = F)) %>%
  gather(., "admix_pop", "ZAnc", c("mex", "maize"))
  
ref_depth2 <- data.frame(chr = pos$chr,
                         pos = pos$pos,
                         mex = apply(mex_ind_anc, 2, mean),
             maize = apply(maize_ind_anc, 2, mean)) %>%
  gather(., "admix_pop", "freq", c("mex", "maize")) %>%
  left_join(ref_depth, by = c("chr", "pos", "admix_pop"))


# plot all barplots together
p <- ggplot(ref_depth, aes(y = ZAnc, x = tot_depth_bin,
                      fill = admix_pop)) + #, color = admix_pop)) +#, color = ZAnc_group)) +
  #scale_colour_manual(values = colors_maize2mex[c(1,4)]) +
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_boxplot() +
  labs(x = "depth of coverage for reference panel", 
       y = "ZAnc statistic by group",
       main = "Effect of maize and mexicana ref. panel depth of coverage on Zanc statistic") +
  facet_wrap(~ ref) + 
  ggtitle("Maize (L) and Mex (R) reference panel depth effect on ZAnc in maize (blue) and in mex. (yellow) admixed pops")
plot(p)
# save plot
ggsave("../plots/ZAnc_in_maize_and_mex_by_ref_panel_depth_cov.png", plot = p, device = png(), 
       width = 20, height = 10, units = "in",
       dpi = 200)

# plot for mean (across ind's not pops) ancestry frequency
# against depth of coverage for reference pops (instead of ZAnc)
p2 <- ggplot(ref_depth2, aes(y = freq, x = tot_depth_bin,
                           fill = admix_pop)) + ##, color = admix_pop)) +
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  #scale_colour_manual(values = colors_maize2mex[c(1,4)]) +
  geom_boxplot() +
  labs(x = "depth of coverage for reference panel", 
       y = "Mex. ancestry freq. by subgroup") +
  facet_wrap(~ ref) + 
  ggtitle("Maize (L) and Mex (R) reference panel depth effect on frequency of mexicana ancestry in maize (blue) and in mex. (yellow) admix pops")
plot(p2)
# save plot
ggsave("../plots/anc_freq_in_maize_and_mex_by_ref_panel_depth_cov.png", 
       plot = p2, device = png(), 
       width = 20, height = 10, units = "in",
       dpi = 200)


# plot one highly significant region and individual ancestry & ref coverage

d1 %>%
  filter(., chr == "4") %>%
  #filter(., pos > 23000000 & pos < 28000000) %>%
  filter(., pos > 6000000 & pos < 10000000) %>%
  # calculate rolling (or windowed) mean anc freq
  mutate(., ZAnc_maize_wind = zoo::rollmean(ZAnc_maize, w, fill = list(NA, NULL, NA))) %>%
  .[c(T, rep(F, as.integer(w/5) - 1)), ] %>% # plot every 5th point
  ggplot(., aes(x = pos, y = ZAnc_maize)) + 
  geom_line(size = .5, alpha = .5, col = "orange")
d1 %>%
  filter(., chr == "4") %>%
  #filter(., pos > 23000000 & pos < 28000000) %>%
  filter(., pos > 6000000 & pos < 10000000) %>%
  # calculate rolling (or windowed) mean anc freq
  mutate(., ZAnc_maize_wind = zoo::rollmean(ZAnc_maize, w, fill = list(NA, NULL, NA))) %>%
  .[c(T, rep(F, as.integer(w/5) - 1)), ] %>% # plot every 5th point
  ggplot(., aes(x = pos, y = ZAnc_maize)) + 
  geom_line(size = .5, alpha = .5, col = "darkred")
region <- which(d1$pos > 6000000 & d1$pos < 10000000 & d1$chr=="4")
plot(y = maize_ind_anc[1, region], 
     x = d1$pos[region], type = "l",
     xlab = "position chr4",
     ylab = "ind. ancestry")
for (i in 1:nrow(maize_ind_anc)){
  lines(y = maize_ind_anc[i, region], 
        x = d1$pos[region])
}
for (i in 1:5){
  lines(y = maize_ind_anc[i, region], 
        x = d1$pos[region])
}
str(maize_ind_anc)

#abline(loess(y = input366$mex_tot[c(T, rep(F, 1000))] ~ maize$ZAnc[c(T, rep(F, 1000))]), 
#col = "blue")
# Look at a few outlier regions in maize --
# what do individual ancestries look like?
# what do posteriors look like?
# what does coverage look like (from ancestry_hmm input file)?
# what is snp density like?
# what does REFERENCE panel coverage look like (from ancestry_hmm input)?

# Now repeat these sanity checks with a few outlier regions in mexicana --

# are there associations between ZAnc and environment?
# start w/ maize

# how does this compare to just plain Anc ~ Env?





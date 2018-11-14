library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(zoo) # for rolling mean (mean across windows)
# helper file with some useful functions
source("../../covAncestry/forqs_sim/k_matrix.R") # import useful functions

#dir_in = "../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/"
dir_in = "../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/"
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
meta <- read.table("../data/pass1_ids.txt", stringsAsFactors = F, 
                   header = T, sep = "\t") %>%
  filter(., est_coverage >= 0.05) %>%
  left_join(., alphas, by = c("popN"))
# individual by individual K matrix
included_inds = meta %>%
  filter(alpha_maize > 0 & alpha_mex > 0) %>%
  .[order(.$popN), ]
pops = unique(included_inds[order(included_inds$popN), c("popN", "zea", "LOCALITY", "alpha_maize", "alpha_mex")])
#write.table(pops$popN, paste0(dir_in, "/included_pops.txt"),
#            col.names = F, row.names = F, quote = F)

# read in population ancestry input files from calc_genomewide_pop_anc_freq.R
pop_anc_list = lapply(pops$popN, function(pop) read.table(paste0(dir_anc, "/pop", pop, ".anc.freq"), 
                                                stringsAsFactors = F))

# combine population ancestry frequencies into a matrix
# where rows are populations and columns are snps
all_anc = t(do.call(cbind, pop_anc_list))
rownames(all_anc) <- paste(substr(pops$zea, 1, 4), substr(pops$LOCALITY, 1, 2), sep = ".")
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
        RowSideColors = col_pop_location)
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
        RowSideColors = col_ind_location)
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
pos <- read.table(
            paste0(dir_in, "input/pos_recomb_rates.txt"), 
            header = T, stringsAsFactors = F)
# use if no recombination rates have been calculated
#pos = read.table(paste0(dir_in, "/output_noBoot/HILO1.posterior"),
#                 stringsAsFactors = F, header = T)[ , 1:2] %>%
#  rename(pos = position)

hist(pos$rate)
summary(pos$rate)
# quantiles for defining 'low' and 'high' recomb.
low_quant = .1
high_quant = .9


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
        RowSideColors = col_pop_location)
dev.off()

png(paste("../plots/K_covplot_pop_combined_high_recomb_rate.png"), # saves plot as pin in ../plots/
    height = 7, width = 7, units = "in", res = 150)
heatmap(high_r$K, Colv = NA, Rowv = NA, main = "high recomb - pop anc cov.",
        col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
        RowSideColors = col_pop_location)
dev.off()


# low vs. high recomb. ind. level K
png(paste("../plots/K_covplot_ind_combined_low_recomb_rate.png"), # saves plot as pin in ../plots/
    height = 14, width = 14, units = "in", res = 150)
heatmap(low_r_ind$K, Colv = NA, Rowv = NA, main = "low recomb - ind anc cov.",
        col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
        RowSideColors = col_ind_location)
dev.off()

png(paste("../plots/K_covplot_ind_combined_high_recomb_rate.png"), # saves plot as pin in ../plots/
    height = 14, width = 14, units = "in", res = 150)
heatmap(high_r_ind$K, Colv = NA, Rowv = NA, main = "high recomb - ind anc cov.",
        col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
        RowSideColors = col_ind_location)
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
# I am not sure whether this is the appropriate
# alpha to recalculate it for just high or low recombining regions --
# I think I should be using a global alpha from genomewide
# or maybe weighing it differently than just taking mean across thinned snps
# because SNPs kept are more likely in low recomb. areas

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
sum(abs(maize[["ZAnc"]]) < maize_bound)/ncol(maize_anc) # 3.6% genome overall appears to be under selection
sum(maize[["ZAnc"]] > maize_bound)/ncol(maize_anc) # nearly all of it in the teosinte direction
mex_bound = qnorm(p = .01, mean = 0, sd = sqrt(nrow(mex_anc)), lower.tail = F)
all_bound = qnorm(p = .01, mean = 0, sd = sqrt(nrow(all_anc)), lower.tail = F)
sum(mex[["ZAnc"]] > mex_bound)/ncol(mex_anc) # hitting boundary?
sum(mex[["ZAnc"]] < -mex_bound)/ncol(mex_anc) # 1-2% maize biased


# histogram of ZAnc statistic
pop_results = list(all, maize, mex)
for (i in pop_results){
  png(paste("../plots/ZAnc_hist_", deparse(substitute(i)) , ".png"), # saves plot as pin in ../plots/
      height = 5, width = 5, units = "in", res = 150)
  hist(i[["ZAnc"]], main = paste("ZAnc", deparse(substitute(i))))
  mtext(paste0("p<.01 outliers maize:", round(sum(i[["ZAnc"]] < -all_bound)/ncol(all_anc), 3), 
               " mex:", round(sum(i[["ZAnc"]] > all_bound)/ncol(all_anc), 3))) # I suspect outliers are driven mainly by maize
  abline(v = c(all_bound, -all_bound), col = "red")
  dev.off()
}


png(paste("../plots/ZAnc_hist_combined.png"), # saves plot as pin in ../plots/
    height = 5, width = 5, units = "in", res = 150)
hist(all[["ZAnc"]], main = "Combined ZAnc")
mtext(paste0("p<.01 outliers maize:", round(sum(all[["ZAnc"]] < -all_bound)/ncol(all_anc), 3), 
                  " mex:", round(sum(all[["ZAnc"]] > all_bound)/ncol(all_anc), 3))) # I suspect outliers are driven mainly by maize
abline(v = c(all_bound, -all_bound), col = "red")
dev.off()
png(paste("../plots/ZAnc_hist_maize.png"), # saves plot as pin in ../plots/
    height = 5, width = 5, units = "in", res = 150)
hist(maize[["ZAnc"]], main = "maize ZAnc")
mtext(paste0("p<.01 outliers maize:", round(sum(maize[["ZAnc"]] < -maize_bound)/ncol(maize_anc), 3), 
             " mex:", round(sum(maize[["ZAnc"]] > maize_bound)/ncol(maize_anc), 3)))
abline(v = c(maize_bound, -maize_bound), col = "red")
dev.off()
png(paste("../plots/ZAnc_hist_mex.png"), # saves plot as pin in ../plots/
    height = 5, width = 5, units = "in", res = 150)
hist(mex[["ZAnc"]], main = "mex ZAnc")
mtext(paste0("p<.01 outliers maize:", round(sum(mex[["ZAnc"]] < -mex_bound)/ncol(mex_anc), 3), 
             " mex:", round(sum(mex[["ZAnc"]] > mex_bound)/ncol(mex_anc), 3)))
abline(v = c(mex_bound, -mex_bound), col = "red") # may be too close to the boundary to get high mex. outliers
dev.off()
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
  geom_point(cex= .1, alpha = .5, col = "orange") + facet_wrap(~chr)
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



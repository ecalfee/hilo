library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(zoo) # for rolling mean (mean across windows)
# helper file with some useful functions
source("../../covAncestry/forqs_sim/k_matrix.R") # import useful functions

dir_in = "../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/"
dir_anc = paste0(dir_in, "output_noBoot/anc/")
# color palette for locations (n=13 locations)
location_colors = c("gray10", "deepskyblue1",
  brewer.pal(11, "Spectral"))
# admixture proportions per pop
alphas <- read.table("../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/input/globalAdmixtureByPopN.txt",
                     stringsAsFactors = F, header = F)
colnames(alphas) <- c("popN", "alpha_maize", "alpha_mex")
# metadata for each individual, including pop association
meta <- read.table("../data/pass1_ids.txt", stringsAsFactors = F, 
                   header = T, sep = "\t") %>%
  filter(est_coverage >= 0.05) %>%
  left_join(., alphas, by = c("popN"))
# individual by individual K matrix
included_inds = meta %>%
  filter(alpha_maize > 0 & alpha_mex > 0) %>%
  .[order(.$popN), ]
pops = unique(included_inds[order(included_inds$popN), c("popN", "zea", "LOCALITY", "alpha_maize", "alpha_mex")])
write.table(pops$popN, paste0(dir_in, "/included_pops.txt"),
            col.names = F, row.names = F, quote = F)

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
pos_thin <- read.table(pos_thin, 
            paste0(dir_in, "input/pos_recomb_rates.txt"), 
            header = T, stringsAsFactors = F)
hist(pos_thin$rate)
summary(pos_thin$rate)
# quantiles for defining 'low' and 'high' recomb.
low_quant = .1
high_quant = .9

# low recombination rate regions <= low quantile:
low_r_anc = all_anc[ , pos_thin$rate <= quantile(pos_thin$rate, low_quant)]
low_r = make_calcs(low_r_anc)
low_r_ind_anc = all_ind_anc[ , pos_thin$rate <= quantile(pos_thin$rate, low_quant)]
low_r_ind = make_calcs(low_r_ind_anc)
# high recombination rate regions >= high quantile:
high_r_anc = all_anc[ , pos_thin$rate >= quantile(pos_thin$rate, high_quant)]
high_r = make_calcs(high_r_anc)
high_r_ind_anc = all_ind_anc[ , pos_thin$rate >= quantile(pos_thin$rate, high_quant)]
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

# Sanity checks:
# Look at a few outlier regions in maize --
# what do individual ancestries look like?
# what do posteriors look like?
# what does coverage look like (from ancestry_hmm input file)?
# what is snp density like?
# what does REFERENCE panel coverage look like (from ancestry_hmm input)?

# Now repeat these sanity checks with a few outlier regions in mexicana --


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

# now get positions so I can plot ZAnc by chromosome:
pos = read.table(paste0(dir_in, "/output_noBoot/HILO1.posterior"),
                   stringsAsFactors = F, header = T)[ , 1:2]
d = data.frame(pos, ZAnc_maize = maize[["ZAnc"]], ZAnc_mex = mex[["ZAnc"]], ZAnc_all = all[["ZAnc"]])
d %>%
  #filter(., chrom == "4") %>%
  ggplot(., aes(x = position, y = ZAnc_maize)) + 
  geom_point(cex= .1, alpha = .5, col = "orange") + facet_wrap(~chrom)
ggsave("../plots/ZAnc_maize_whole_genome.png", height = 15, width = 15, units = "in", device = "png")
d %>%
  #filter(., chrom == "4") %>%
  ggplot(., aes(x = position, y = ZAnc_mex)) + 
  geom_point(cex= .1, alpha = .5, col = "blue") + facet_wrap(~chrom)
ggsave("../plots/ZAnc_mex_whole_genome.png", height = 15, width = 15, units = "in", device = "png")
d %>%
  #filter(., chrom == "4") %>%
  ggplot(., aes(x = position, y = ZAnc_all)) + 
  geom_point(cex= .1, alpha = .5, col = "black") + facet_wrap(~chrom)
ggsave("../plots/ZAnc_combined_whole_genome.png", height = 15, width = 15, units = "in", device = "png")


# Now get smoothed ZAnc statistic plots 
# (taking mean across windows of 50 thinned positons ~ 50-100 kb)
# geom_smooth (method = "gam" (loess is too memory intense for > 1000 points))
# geom_smooth uses a regression method (of your choice) ..
# to actually just use raw mean of a window I'll use roll_mean() from package zoo
w = 100 # window size (# thinned snps)
d %>%
  #filter(., chrom == "4") %>%
  # calculate rolling (or windowed) mean ZAnc
  mutate(., ZAnc_maize_wind = zoo::rollmean(ZAnc_maize, w, fill = list(NA, NULL, NA))) %>%
  .[c(T, rep(F, as.integer(w/5) - 1)), ] %>% # plot every 5th point
  ggplot(., aes(x = position, y = ZAnc_maize)) + 
  geom_line(size = .5, cex= .1, alpha = .5, col = "orange") + 
  #geom_line(size = 1.2, cex= .1, alpha = .5, col = "orange") +
  facet_wrap(~chrom)

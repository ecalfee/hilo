# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(reshape2)
library(VennDiagram)
library(ggrastr)
library(ggupset) # to plot all combinations of outlier sharing across pops

fdr_range = c(0.1, .05, .01)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# TO DO:
# write out generalized least squares comparison to what I'm doing. IF I use AIC it should be AICc for n=14 data points.

rerun_all_models <- F # rerun or just load results of models

# load external scripts
source("k_matrix.R") # import useful functions

# load data 
PREFIX="pass2_alloMAIZE"

# metadata for each individual, including pop association
pop_elev <- read.table("../data/riplasm/gps_and_elevation_for_sample_sites.txt",
                       stringsAsFactors = F, header = T, sep = "\t") %>%
  dplyr::select(., popN, ELEVATION, LAT, LON)
hilo <- read.table("../samples/hilo_meta.txt", stringsAsFactors = F, header = T, sep = "\t")
# admixture proportions per pop
meta.pops <- read.table(paste0("../global_ancestry/results/NGSAdmix/", PREFIX, "/globalAdmixtureIncludedByPopN.txt"),
                        stringsAsFactors = F, header = T) %>%
  left_join(., pop_elev, by = c("popN")) %>%
  left_join(., unique(dplyr::select(hilo, c("popN", "zea", "symp_allo", "RI_ACCESSION", "GEOCTY", "LOCALITY"))), by = c("popN")) %>%
  mutate(pop = paste0("pop", popN))

# genomic sites information
# SNPs with ancestry calls
dir_sites = paste0("../local_ancestry/results/thinnedSNPs/", PREFIX)
sites <- do.call(rbind, 
                 lapply(1:10, function(i)
                   read.table(paste0(dir_sites, "/chr", i, ".var.sites"),
                              header = F, stringsAsFactors = F)))
colnames(sites) <- c("chr", "pos", "major", "minor")
# get inversions:
inv = read.table("../data/refMaize/inversions/knownInv_v4_coord.txt",
                 stringsAsFactors = F, header = T)
colnames(inv) = c("ID", "chr", "start", "end", "length")
excl_inv_buffer = 1000000 # exclude 1 megabase around inversion
inv$excl_start = inv$start - excl_inv_buffer
inv$excl_end = inv$end + excl_inv_buffer
# which sites are within the inversion?
sites <- mutate(sites, inv4m_1Mb_buffer = (chr == inv$chr[inv$ID=="inv4m"] & 
                                             pos >= inv$excl_start[inv$ID=="inv4m"] & 
                                             pos <= inv$excl_end[inv$ID == "inv4m"]),
                inv4m = (chr == inv$chr[inv$ID=="inv4m"] & 
                           pos >= inv$start[inv$ID=="inv4m"] & 
                           pos <= inv$end[inv$ID == "inv4m"])) 


# ancestry calls
LOCAL_ANC_SUBDIR="output_noBoot"
dir_in = paste0("../local_ancestry/results/ancestry_hmm/", PREFIX)
dir_anc = file.path(dir_in, LOCAL_ANC_SUBDIR, "anc")
# read in population ancestry input files from calc_genomewide_pop_anc_freq.R
pop_anc_list = lapply(meta.pops$popN, function(pop) read.table(paste0(dir_anc, "/pop", pop, ".anc.freq"), 
                                                               stringsAsFactors = F))

# combine population ancestry frequencies into a matrix
# where rows are populations and columns are snps
all_anc = do.call(cbind, pop_anc_list)
colnames(all_anc) <- meta.pops$pop
rm(pop_anc_list)

# get mean ancestry across individuals sampled (NOT mean of pops)
maize_pops <- unique(meta.pops$pop[meta.pops$zea == "maize"])
mexicana_pops <- unique(meta.pops$pop[meta.pops$zea == "mexicana"])
maize_anc <- all_anc[ , maize_pops]
mexicana_anc <- all_anc[ , mexicana_pops]
anc <- sites
anc$pop_meanAlpha_maize = rowMeans(all_anc[ , maize_pops])
anc$pop_meanAlpha_mex = rowMeans(all_anc[ , mexicana_pops])
N_per_pop_mex <- sapply(mexicana_pops, function(pop) unique(meta.pops$n[meta.pops$pop == pop]))
anc$meanAlpha_mex = t((N_per_pop_mex %*% t(all_anc[ , mexicana_pops]))/sum(N_per_pop_mex))[ , 1]
N_per_pop_maize <- sapply(mexicana_pops, function(pop) unique(meta.pops$n[meta.pops$pop == pop]))
anc$meanAlpha_maize = t((N_per_pop_maize %*% t(all_anc[ , maize_pops]))/sum(N_per_pop_maize))[ , 1]
# file for Jeff
file4jeff <- left_join(anc, d_simple_bElev,
                       by = colnames(d_simple_bElev)[1:6]) %>%
  dplyr::select(c(colnames(d_simple_bElev), "meanAlpha_maize", "meanAlpha_mex")) %>%
  rename(slope_elev_sig = significance) %>%
  mutate(meanAlpha_maize = round(meanAlpha_maize, 3),
         meanAlpha_mex = round(meanAlpha_mex, 3),
         slope_elev = round(slope_elev, 6))
write.table(file4jeff, "results/outliers_for_jeff_9.5.19.txt",
            col.names = T, row.names = F, quote = F, sep= "\t")

# calculations
zAnc_mexicana = make_K_calcs(t(mexicana_anc))
zAnc_maize = make_K_calcs(t(maize_anc))

write.table(zAnc_maize$K, "results/zAnc_K_matrix_maize.txt", sep = "\t", row.names = T, col.names = T, quote = F)
write.table(zAnc_maize$alpha, "results/zAnc_alpha_maize.txt", sep = "\t", row.names = T, col.names = T, quote = F)

# script to run test cases of ancestry selection models using empirical K matrices
#source("zanc_selection_models.R") # makes plots

# calculate zAnc at each locus. 
zAnc3 <- t(apply(maize_anc, 1, function(l) zAnc(ancFreq = l, 
                                                invL = zAnc_maize$InvL, 
                                                alpha = zAnc_maize$alpha)))
colnames(zAnc3) <- colnames(maize_anc)

zAnc3_mex <- t(apply(mexicana_anc, 1, function(l) zAnc(ancFreq = l, 
                                                       invL = zAnc_mexicana$InvL, 
                                                       alpha = zAnc_mexicana$alpha)))
colnames(zAnc3_mex) <- colnames(mexicana_anc)
# Null model
# Now get zTz -- basically the sum of residual errors to the null model
# zTz
zTz3 <- apply(zAnc3, 1, function(i) t(i) %*% i)
zTz3_mex <- apply(zAnc3_mex, 1, function(i) t(i) %*% i)


# Environmental selection models (in maize):
# Elevation-based selection
if (rerun_all_models){ # rerun or just load results of models
  zBeta3_elev_maize <- apply(maize_anc, 1, function(l) # calculate slope for each locus
    zBeta3(ancFreq = l,
           # center elevation
           envWeights = meta.pops$ELEVATION[meta.pops$zea == "maize"] - mean(meta.pops$ELEVATION[meta.pops$zea == "maize"]), 
           invL = zAnc_maize$InvL, 
           alpha = zAnc_maize$alpha,
           zInt = T)) # has transformed intercept
  zb3_elev <- data.frame(t(zBeta3_elev_maize)) # data frame is easier to work with
  write.table(zb3_elev,
              paste0("results/models/", PREFIX, "/maize/elevation.txt"),
              sep = "\t",
              quote = F,
              col.names = T, row.names = F)
}else{
  zb3_elev <- bind_cols(sites, read.table(paste0("results/models/", PREFIX, "/maize/elevation.txt"),
                         sep = "\t",
                         header = T,
                         stringsAsFactors = F))
}

# 
zb3_elev_noncentered <- bind_cols(sites, read.table(paste0("results/models/", PREFIX, "/maize/elevation_noncentered.txt"),
                                        sep = "\t",
                                        header = T,
                                        stringsAsFactors = F))

# threshold Elevation-based selection - 11pops > 1900m experience (equal) selection, 3 lowland pops do not
if (rerun_all_models){ # rerun or just load results of models
  zBeta3_over1900m_maize <- apply(maize_anc, 1, function(l) # calculate slope for each locus
    zBeta3(ancFreq = l,
           envWeights = as.numeric(meta.pops$ELEVATION[meta.pops$zea == "maize"] > 1900), 
           invL = zAnc_maize$InvL, 
           alpha = zAnc_maize$alpha,
           zInt = F)) # no intercept -- low elevation populations not under selection
  zb3_over1900m <- data.frame(t(zBeta3_over1900m_maize)) # data frame is easier to work with
  write.table(zb3_over1900m,
              paste0("results/models/", PREFIX, "/maize/sel_over1900m.txt"),
              sep = "\t",
              quote = F,
              col.names = T, row.names = F)
}else{
  zb3_over1900m <- bind_cols(sites, read.table(paste0("results/models/", PREFIX, "/maize/sel_over1900m.txt"),
                         sep = "\t",
                         header = T,
                         stringsAsFactors = F))
}
# is elevation estimated positive or negative for top outliers?
zb3_elev %>%
  bind_cols(., sites) %>%
  mutate(zTz3 = zTz3) %>%
  ggplot(aes(x = log10(pval_zEnv), y = log10(zTz3), color = zEnv)) +
  geom_point()
# take subset to plot
# effect of centering environment should be to just change interpretation of intercept
# as effect of elevation on mexicana ancestry at the mean of our sampled elevations
# slope with environment doesn't change
plot(zb3_elev_noncentered$zEnv[c(T, rep(F, 100))], zb3_elev$zEnv[c(T, rep(F, 100))])
plot(zb3_elev_noncentered$pval_zEnv[c(T, rep(F, 100))], zb3_elev$pval_zEnv[c(T, rep(F, 100))])
# intercept of course changes
plot(zb3_elev_noncentered$zInt[c(T, rep(F, 100))], zb3_elev$zInt[c(T, rep(F, 100))])
plot(zb3_elev_noncentered$pval_zInt[c(T, rep(F, 100))], zb3_elev$pval_zInt[c(T, rep(F, 100))])
# centering, as expected, gets rid of strong negative relationship 
# between estimated intercept and slope with environment
plot(zb3_elev_noncentered$zInt[c(T, rep(F, 100))], zb3_elev_noncentered$zEnv[c(T, rep(F, 100))])
plot(zb3_elev$zInt[c(T, rep(F, 100))], zb3_elev$zEnv[c(T, rep(F, 100))])
# intercept may be more predictive than ~ elevation for a general misfit to the neutral model
# in particular, high positive intercepts don't fit well
plot(zTz3[c(T, rep(F, 100))], zb3_elev$zInt[c(T, rep(F, 100))])
plot(zTz3[c(T, rep(F, 100))], zb3_elev$zEnv[c(T, rep(F, 100))])
contour(x = zb3_elev$zInt[c(T, rep(F, 100))], # would need to be sorted low to high
        y = zb3_elev$zEnv[c(T, rep(F, 100))],
        z = zTz3[c(T, rep(F, 100))]) # doesn't work in r to visualize this way

hist(zTz3)
# ~zElev individual loci
bind_cols(maize_anc, sites) %>%
  bind_cols(., zb3_elev) %>%
  mutate(zTz3 = zTz3) %>%
  mutate(elev_pval = ifelse(pval_zEnv < 0.01, "zElev pval<.01", "pval>=.01")) %>%
  ggplot(aes(x = zTz3, fill = zEnv>0)) +
  geom_histogram(position = "stack") +
  facet_wrap(~elev_pval, scales = "free_y")+
  ggtitle("zTz scores for low vs. high p-val from zAnc~zElev slope")
ggsave("plots/zEnv_pval_elevation_vs_ztz_score_all_loci.png", device = "png",
       width = 10, height = 6, units = "in")
# plot the inverse
bind_cols(maize_anc, sites) %>%
  bind_cols(., zb3_elev) %>%
  mutate(zTz3 = zTz3) %>%
  mutate(elev_pval = ifelse(zTz3 > quantile(zTz3, .99), "zTz > 99%", "zTz non-outlier")) %>%
  ggplot(aes(x = log10(pval_zEnv), fill = zEnv>0)) +
  geom_histogram(position = "stack") +
  facet_wrap(~elev_pval, scales = "free_y")+
  ggtitle("p-vals slopes from zAnc~zElev for low vs. high zTz loci")
ggsave("plots/zEnv_pval_elevation_vs_99percent_ztz_score_all_loci.png", device = "png",
       width = 10, height = 6, units = "in")


# make plot log 10 pval for zEnv vs. zElev slope estimate to show enrichment for selection for mexicana ancestry at high elevation
zb3_elev %>%
  bind_cols(., sites) %>%
  ggplot(aes(x = zEnv, y = -log10(pval_zEnv), color = inv4m)) +
  geom_point(alpha = .2) +
  ggtitle("model: zAnc ~ elevation")
ggsave("plots/zbe_elev_pval_vs_effect.png", device = "png",
       height = 6, width = 8, units = "in")
png("plots/zb3_elev_pval_combined.png")
hist(zb3_elev$pval_zEnv)
dev.off()
hist(zb3_elev$pval_zEnv[zb3_elev$zEnv >= 0])
hist(zb3_elev$pval_zEnv[zb3_elev$zEnv <= 0])
ggplot(zb3_elev, aes(x = pval_zEnv, fill = (zEnv >= 0))) + 
  geom_histogram(binwidth = .01) +
  facet_wrap(~(zEnv >= 0)) +
  ggtitle("pvalue distribution zAnc ~ zElev, TRUE = + slope")
ggsave("plots/zbe_elev_pval_histograms.png", device = "png",
       height = 6, width = 8, units = "in")

# Universal selection for or against mexicana
if (rerun_all_models){ # rerun or just load results of models
  zBeta3_highmex_maize <- apply(maize_anc, 1, function(l) # calculate slope for each locus
    zBeta3(ancFreq = l,
           envWeights = rep(1, 14), 
           invL = zAnc_maize$InvL, 
           alpha = zAnc_maize$alpha,
           zInt = F)) # no intercept 
  zb3_hmex <- data.frame(t(zBeta3_highmex_maize)) # data frame is easier to work with
  write.table(zb3_hmex,
              paste0("results/models/", PREFIX, "/maize/universal_sel_hmex.txt"),
              sep = "\t",
              quote = F,
              col.names = T, row.names = F)
}else{
  zb3_hmex <- bind_cols(sites, read.table(paste0("results/models/", PREFIX, "/maize/universal_sel_hmex.txt"),
                         sep = "\t",
                         header = T,
                         stringsAsFactors = F))
}
# compare intercept in elevation model to slope on 'all pops' model. 
# Why aren't these the same? Is this just fitting error?
plot(zb3_elev$zInt[c(T, rep(F, 100))], zb3_hmex$zEnv[c(T, rep(F, 100))])
plot(zb3_elev$pval_zInt[c(T, rep(F, 100))], zb3_hmex$pval_zEnv[c(T, rep(F, 100))])
# should be strongly negatively correlated
plot(zb3_elev$sum_sq_res[c(T, rep(F, 100))] - zb3_hmex$sum_sq_res[c(T, rep(F, 100))],
     zb3_elev$pval_zEnv[c(T, rep(F, 100))],
     col = ifelse(sites$inv4m[c(T, rep(F, 100))], "blue", "black")) # mostly makes sense
plot(1 - (zb3_elev$sum_sq_res[c(T, rep(F, 100))]/zb3_hmex$sum_sq_res[c(T, rep(F, 100))]),
     zb3_elev$pval_zEnv[c(T, rep(F, 100))],
     col = ifelse(sites$inv4m[c(T, rep(F, 100))], "blue", "black")) # mostly makes sense
plot(zb3_hmex$sum_sq_res[c(T, rep(F, 100))] - zTz3[c(T, rep(F, 100))],
     zb3_hmex$pval_zEnv[c(T, rep(F, 100))],
     col = ifelse(sites$inv4m[c(T, rep(F, 100))], "blue", "black"))
plot(1 - (zb3_hmex$sum_sq_res[c(T, rep(F, 100))]/zTz3[c(T, rep(F, 100))]),
     zb3_hmex$pval_zEnv[c(T, rep(F, 100))],
     col = ifelse(sites$inv4m[c(T, rep(F, 100))], "blue", "black"))
plot(zb3_elev$sum_sq_res[c(T, rep(F, 100))] - zTz3[c(T, rep(F, 100))],
     zb3_elev$zEnv[c(T, rep(F, 100))],
     col = ifelse(sites$inv4m[c(T, rep(F, 100))], "blue", "black"))
plot(zb3_elev$sum_sq_res[c(T, rep(F, 100))]-zTz3[c(T, rep(F, 100))],
     zb3_elev$zInt[c(T, rep(F, 100))],
     col = ifelse(sites$inv4m[c(T, rep(F, 100))], "blue", "black"))
plot(1 - (zb3_elev$sum_sq_res[c(T, rep(F, 100))]/zTz3[c(T, rep(F, 100))]),
     zb3_elev$pval_zEnv[c(T, rep(F, 100))],
     col = ifelse(sites$inv4m[c(T, rep(F, 100))], "blue", "black"))
plot(1 - (zb3_elev$sum_sq_res[c(T, rep(F, 100))]/zTz3[c(T, rep(F, 100))]),
     zb3_elev$pval_zInt[c(T, rep(F, 100))],
     col = ifelse(sites$inv4m[c(T, rep(F, 100))], "blue", "black"))
quantile(zb3_elev$sum_sq_res - zTz3, .01)

png("plots/zb3_hmex_pval_combined.png")
hist(zb3_hmex$pval_zEnv)
dev.off()
hist(zb3_hmex$pval_zEnv[zb3_hmex$zEnv >= 0])
hist(zb3_hmex$pval_zEnv[zb3_hmex$zEnv <= 0])
ggplot(zb3_hmex, aes(x = pval_zEnv, fill = (zEnv >= 0))) + 
  geom_histogram(binwidth = .01) +
  facet_wrap(~(zEnv >= 0)) +
  ggtitle("pvalue distribution zAnc ~ z1's, all selected, TRUE = + slope")
ggsave("plots/zbe_hmex_pval_histograms.png", device = "png",
       height = 6, width = 8, units = "in")
zb3_hmex %>%
  bind_cols(., sites) %>%
  ggplot(aes(x = zEnv, y = -log10(pval_zEnv), color = inv4m)) +
  geom_point(size = .1) +
  ggtitle("model: zAnc ~ 1")
ggsave("plots/zbe_hmex_pval_vs_effect.png", device = "png",
       height = 6, width = 8, units = "in")


# test
zb3_onepop_test <- lapply(1:2, function(x)
  data.frame(t(apply(maize_anc[1:10,], 1, function(l) # calculate slope for each locus
    zBeta3(ancFreq = l,
           envWeights = onePop[x, ], # column for 1 population with selection 
           invL = zAnc_maize$InvL, 
           alpha = zAnc_maize$alpha,
           zInt = F))), stringsAsFactors = F)) # no intercept


# Selection in one population only (14 possibilities for populations)
onePop <- matrix(0, 14, 14)
diag(onePop) <- 1 # matrix with all posibilities for which population has selection '1'
if (rerun_all_models){ # rerun or just load results of models
  for (i in 1:14){
    zb3_onepop1 <- apply(maize_anc, 1, function(l) # calculate slope for each locus
        zBeta3(ancFreq = l,
               envWeights = onePop[i, ], # column for 1 population with selection 
               invL = zAnc_maize$InvL, 
               alpha = zAnc_maize$alpha,
               zInt = F))
    zb3_onepop2 <- data.frame(t(zb3_onepop1), stringsAsFactors = F) # no intercept

  write.table(zb3_onepop2,
              paste0("results/models/", PREFIX, "/maize/sel_pop_", i, "_only.txt"),
              sep = "\t",
              quote = F,
              col.names = T, row.names = F)
  }
}
# proportional to the Sum of Squared Errors (SSE), and the log likelihood
# load model results for all 14 models
zb3_onepop <- lapply(1:14, function(i)
  bind_cols(sites, read.table(paste0("results/models/", PREFIX, "/maize/sel_pop_", i, "_only.txt"),
             sep = "\t",
             header = T,
             stringsAsFactors = F)))
ordered_maize_pops <- left_join(data.frame(pop = colnames(maize_anc), stringsAsFactors = F), meta.pops, by = "pop")


png("plots/pos_sel_1pop_pvals_hist.png")
par(mfrow = c(4,4))
lapply(1:14, function(i) hist(zb3_onepop[[i]]$pval_zEnv[zb3_onepop[[i]]$zEnv  > 0], main = paste("+ sel", ordered_maize_pops$LOCALITY[i],
                                                                          ordered_maize_pops$ELEVATION[i], "m"), xlab = "pval", ylab = "freq"))
par(mfrow = c(1,1))
dev.off()

png("plots/neg_sel_1pop_pvals_hist.png")
par(mfrow = c(4,4))
lapply(1:14, function(i) hist(zb3_onepop[[i]]$pval_zEnv[zb3_onepop[[i]]$zEnv  < 0], main = paste("- sel", ordered_maize_pops$LOCALITY[i],
                                                                                                 ordered_maize_pops$ELEVATION[i], "m"), xlab = "pval", ylab = "freq"))
par(mfrow = c(1,1))
dev.off()

with(zb3_elev, plot(r_squared ~ zEnv, col = ifelse(d$inv4m, "blue", "black"), main = "zBeta3 association mexicana ancestry with high elevation in maize"))
with(zb3_elev, plot(sum_sq_res ~ zEnv, col = ifelse(d$inv4m, "blue", "black"), main = "zBeta3 association mexicana ancestry with high elevation in maize"))
with(zb3_hmex, plot(r_squared ~ zEnv, col = ifelse(d$inv4m, "blue", "black"), main = "zBeta3 universally high mexicana anc in maize"))
with(zb3_hmex, plot(sum_sq_res ~ zEnv, col = ifelse(d$inv4m, "blue", "black"), main = "zBeta3 universally high mexicana anc in maize"))

# compare models:
zb3_compare <- data.frame(zTz = zTz3,
                          elev = zb3_elev$sum_sq_res,
                          hmex = zb3_hmex$sum_sq_res,
                          over1900m = zb3_over1900m$sum_sq_res)

# combine all model outputs together into one dataframe
zb3_hmex$model = "hmex"
zb3_elev$model = "elevation"
zb3_over1900m$model = "over1900m"
for (i in 1:14){
  zb3_onepop[[i]]$model <- ordered_maize_pops$LOCALITY[i]
}
zb3_all <- rbind(zb3_hmex, zb3_elev, zb3_over1900m, 
                 do.call(rbind,
                         zb3_onepop))

# note: some of these models are nested: zTz < hmex < elev
# but others are separate model families, e.g. zTz < pop1
# AIC in OLS framework: 
# (sources: https://stats.stackexchange.com/questions/261273/how-can-i-apply-akaike-information-criterion-and-calculate-it-for-linear-regress
# https://en.wikipedia.org/wiki/Akaike_information_criterion)
# AIC = n*ln(SSE/n) + 2k where k is predictors + intercept and n is # data points
# AICc = AIC + (2k^2 + 2k)/(n-k-1) correction for small sample size adds additoinal penalty for model complexity
aic <- function(SSE, n, k, fixed_var = F){
    -2*loglik_OLS(SSE, n, fixed_var) + 2*k # general equation
    # n*log(SSE/n) + 2*k # true in the case of not-fixed variance
}
aic_c <- function(SSE, n, k, fixed_var = F){
  aic(SSE, n, k, fixed_var) + (2*k^2 + 2*k)/(n - k - 1)
}
# log likelihood for OLS regression: https://en.wikipedia.org/wiki/Akaike_information_criterion
loglik_OLS <- function(SSE, n, fixed_var = F){
  # normal LL = -n/2*ln(2*pi*var) - 1/(2*var)*SSE
  if (fixed_var) {
    - n/2*log(2*pi) - 1/2*SSE # under fixed variance of the residuals, var = 1
    } else{
      - n/2*log(2*pi) - n/2*log(SSE/n) - n/2 # under point estimation of variance of the residulas, var = SSE/n, which I've plugged in here
    }  
} 


# Akaike weights for model averaging: calc relative likelihood of each model exp(-0.5 * ∆AIC) and divide by sum of these values across all models. 
aic_w <- function(AICs){
  w = exp(-0.5*(AICs - min(AICs)))
  w/sum(w)
  }


# regular AIC
aic_onepop <- do.call(cbind,
             lapply(zb3_onepop, function(x)
               aic(SSE = x$sum_sq_res, n = 14, k = 1)))
colnames(aic_onepop) <- ordered_maize_pops$LOCALITY
aic_compare <- data.frame(zTz = aic(SSE = zb3_compare$zTz, n = 14, k = 0),
                          elev = aic(SSE = zb3_compare$elev, n = 14, k = 2),
                          hmex = aic(SSE = zb3_compare$hmex, n = 14, k = 1),
                          over1900m = aic(SSE = zb3_compare$over1900m, n = 14, k = 1)) %>%
  bind_cols(., data.frame(aic_onepop))
# AICc for small sample size
aic_c_onepop <- do.call(cbind,
                        lapply(zb3_onepop, function(x)
                          aic_c(SSE = x$sum_sq_res, n = 14, k = 1)))
colnames(aic_c_onepop) <- ordered_maize_pops$LOCALITY
aic_c_compare <- data.frame(zTz = aic_c(SSE = zb3_compare$zTz, n = 14, k = 0),
                            elev = aic_c(SSE = zb3_compare$elev, n = 14, k = 2),
                            hmex = aic_c(SSE = zb3_compare$hmex, n = 14, k = 1),
                            over1900m = aic_c(SSE = zb3_compare$over1900m, n = 14, k = 1)) %>%
  bind_cols(., data.frame(aic_c_onepop))

# AICc for small sample size assuming fixed variance = 1 for all residuals
aic_c_onepop_fixedv <- do.call(cbind,
                        lapply(zb3_onepop, function(x)
                          aic_c(SSE = x$sum_sq_res, n = 14, k = 1, fixed_var = T)))
colnames(aic_c_onepop_fixedv) <- ordered_maize_pops$LOCALITY
aic_c_compare_fixedv <- data.frame(zTz = aic_c(SSE = zb3_compare$zTz, n = 14, k = 0, fixed_var = T),
                            elev = aic_c(SSE = zb3_compare$elev, n = 14, k = 2, fixed_var = T),
                            hmex = aic_c(SSE = zb3_compare$hmex, n = 14, k = 1, fixed_var = T),
                            over1900m = aic_c(SSE = zb3_compare$over1900m, n = 14, k = 1, fixed_var = T)) %>%
  bind_cols(., data.frame(aic_c_onepop_fixedv))


# AIC model weights
aic_w_compare <- data.frame(t(apply(aic_compare, 1, aic_w)))
aic_c_w_compare <- data.frame(t(apply(aic_c_compare, 1, aic_w)))
summary(apply(aic_w_compare, 1, max)) # some loci have higher weight to many models vs. 1 model
apply(aic_w_compare, 2, mean) # what are the dominant modes of selection? mean weights of different models across all loci
apply(aic_c_w_compare, 2, mean) # what are the dominant modes of selection? mean weights of different models across all loci
# I don't find this result all that believable
apply(aic_c_w_compare, 2, summary) # null model zTz never wins

# am I calculating this correctly below? bayes factor values are very high
# loglikelihood to compute bayes factor (ratio of likelihood under alt hyp/ likelihood under null hyp)
ll_onepop <- do.call(cbind,
                        lapply(zb3_onepop, function(x)
                          loglik_OLS(SSE = x$sum_sq_res, n = 14)))
colnames(ll_onepop) <- ordered_maize_pops$LOCALITY
ll <- data.frame(zTz = loglik_OLS(SSE = zb3_compare$zTz, n = 14),
                 elev = loglik_OLS(SSE = zb3_compare$elev, n = 14),
                 hmex = loglik_OLS(SSE = zb3_compare$hmex, n = 14),
                 hmex = loglik_OLS(SSE = zb3_compare$over1900m, n = 14)) %>%
  bind_cols(., data.frame(ll_onepop))
# bayes factor is the likelihood(h1)/likelihood(h2),
# so I have to exponentiate after subtracting the log likelihoods of the alt and null hypothesis
bf <- data.frame(apply(ll, 2, function(x) exp(x - ll$zTz)))

table(apply(aic_compare, 1, which.min))/nrow(aic_compare) # not good - MOST loci should fit the null model best
table(apply(aic_compare[ , 1:3], 1, which.min))/nrow(aic_compare) # but if I limit it to just a few models, the neutral model wins most of the time
table(apply(aic_compare[zTz3 > quantile(zTz3, .99),], 1, which.min))/sum(zTz3 > quantile(zTz3, .99))
table(apply(aic_compare[zTz3 > quantile(zTz3, .95),], 1, which.min))/sum(zTz3 > quantile(zTz3, .95))
table(apply(aic_compare[zTz3 > quantile(zTz3, .9),], 1, which.min))/sum(zTz3 > quantile(zTz3, .9))

table(apply(aic_c_compare, 1, which.min))/nrow(aic_c_compare)
tab1 <- table(apply(aic_c_compare, 1, which.min))/nrow(aic_c_compare)
names(tab1) <- colnames(aic_c_compare)
png("plots/counts_model_win_all_loci_18models.png", width = 8, height = 6, res = 300, units = "in")
barplot(tab1, las = 2, main = "all loci - which model has lowest AICc?")
dev.off()

tab2 <- table(apply(aic_c_compare[ , 1:5], 1, which.min))/nrow(aic_c_compare)  # 80% of the time the neutral model does better
names(tab2) <- colnames(aic_c_compare)[1:5]
png("plots/counts_model_win_all_loci_5models.png", width = 8, height = 6, res = 300, units = "in")
barplot(tab2, las = 2, main = "all loci - which model has lowest AICc?")
dev.off()

# fixed variance assumption
tab1_fixedv <- table(apply(aic_c_compare_fixedv, 1, which.min))/nrow(aic_c_compare_fixedv)
names(tab1_fixedv) <- colnames(aic_c_compare_fixedv)
barplot(tab1_fixedv, las = 2, main = "all loci fixed var=1, which model has lowest AICc?")
tab1_fixedv_99 <- table(apply(aic_c_compare_fixedv[zTz3 > quantile(zTz3, .99),], 1, which.min))/nrow(aic_c_compare_fixedv[zTz3 > quantile(zTz3, .99),])
names(tab1_fixedv_99) <- colnames(aic_c_compare_fixedv)[as.numeric(unlist(dimnames(tab1_fixedv_99)))]
barplot(tab1_fixedv_99, las = 2, main = "99% zTz top outlier loci var=1 fixed - which model has lowest AICc?")


# just tabulate for outlier high zTz loci
tab1_99 <- table(apply(aic_c_compare[zTz3 > quantile(zTz3, .99),], 1, which.min))/nrow(aic_c_compare[zTz3 > quantile(zTz3, .99),])
names(tab1_99) <- colnames(aic_c_compare)
png("plots/counts_model_win_99perc_outlier_loci_18models.png", width = 8, height = 6, res = 300, units = "in")
barplot(tab1_99, las = 2, main = "99% zTz top outlier loci - which model has lowest AICc?")
dev.off()

aic_c_winners <- data.frame(x = apply(aic_c_compare, 1, which.min)) %>%
  left_join(., data.frame(top_model = c("zTz", "elevation", "hmex", "over1900m", ordered_maize_pops$LOCALITY), 
                          x = 1:ncol(aic_c_compare), 
                          stringsAsFactors = F), by = "x") %>%
  bind_cols(sites[,c("chr", "pos")], .) %>%
  left_join(., zb3_all, by = c("chr", "pos", "top_model"="model"))
aic_c_winners %>%
  filter(zTz3 > quantile(zTz3, .99)) %>%
  ggplot(aes(x = top_model, fill = zEnv > 0)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("AICc top model (maize); 99% zTz outlier loci only, blue = sel for mexicana")
ggsave("plots/AICc_top_model_99perzTz_outliers.png", device = "png",
       width = 10, height = 6, units = "in")


# (!) it looks like most of the environmental selection is towards low mex at high altitude
# (which makes little sense). But there's also the intercept. Which is not only positive,
# but high for these outlier loci
aic_c_winners %>%
  filter(zTz3 > quantile(zTz3, .99)) %>%
  ggplot(aes(x = top_model, fill = zInt > 1)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("AICc top model (maize); 99% zTz outlier loci only, blue = intercept high towards mexicana")
ggsave("plots/AICc_top_model_99perzTz_outliers_zInt.png", device = "png",
       width = 10, height = 6, units = "in")
# intercept isn't always high:
hist(aic_c_winners$zInt[zTz3 > quantile(zTz3, .99) & aic_c_winners$top_model == "elevation"])
hist(aic_c_winners$zInt[aic_c_winners$top_model == "elevation"])
hist(aic_c_winners$zEnv[zTz3 > quantile(zTz3, .99) & aic_c_winners$top_model == "hmex"])
# so what do these loci look like?
# ~zElev individual loci. The highest zTz3 outliers all have elevated mexicana (high intercepts)
# and flat or mild downward slope. broadening to top 90% percentile and we see some upward 
# trending loci
# loci that are universally selected against mexicana are never in this set of top zTz3 outliers
bind_cols(maize_anc, sites) %>%
  bind_cols(., zb3_elev) %>%
  mutate(zTz3 = zTz3) %>%
  #filter(., aic_c_winners$top_model == "elevation" & zTz3 > quantile(zTz3, .99)) %>%
  filter(., aic_c_winners$top_model == "elevation" & zTz3 > quantile(zTz3, .9)) %>%
  #filter(., aic_c_winners_fixedv$top_model == "elevation" & zTz3 > quantile(zTz3, .9)) %>%
  .[c(T, rep(F, 50)), ] %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = mex_freq, shape = inv4m, 
             color = log10(zTz3), group = pos)) +
  geom_point() +
  geom_line() +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset 90% top zTz3 outliers w/ elevation model best fit")
ggsave(paste0("plots/a_few_example_loci_zElev_best_model_highzTz3_outliers.png"),
       device = "png",
       width = 10, height = 8, units = "in")



aic_c_winners %>%
  ggplot(aes(x = top_model, fill = zEnv > 0)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle("AICc top model (maize); all loci. blue = sel for mexicana")
ggsave("plots/AICc_top_model_all_loci.png", device = "png",
       width = 10, height = 6, units = "in")

aic_c_winners_fixedv <- data.frame(x = apply(aic_c_compare_fixedv, 1, which.min)) %>%
  left_join(., data.frame(top_model = c("zTz", "elevation", "hmex", "over1900m", ordered_maize_pops$LOCALITY), 
                          x = 1:ncol(aic_c_compare_fixedv), 
                          stringsAsFactors = F), by = "x") %>%
  bind_cols(sites[,c("chr", "pos")], .) %>%
  left_join(., zb3_all, by = c("chr", "pos", "top_model"="model")) 

aic_c_winners_fixedv %>%
  filter(zTz3 > quantile(zTz3, .99)) %>%
  ggplot(aes(x = top_model, fill = zEnv > 0)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle("AICc top model (maize); fixed var = 1; 99% zTz outliers, blue = sel for mexicana")
ggsave("plots/AICc_top_model_fixedvar1_99perzTz_outliers.png", device = "png",
       width = 10, height = 6, units = "in")
aic_c_winners_fixedv %>%
  ggplot(aes(x = top_model, fill = zEnv > 0)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle("AICc top model (maize); fixed var = 1; all loci, blue = sel for mexicana")
ggsave("plots/AICc_top_model_fixedvar1_allLoci_outliers.png", device = "png",
       width = 10, height = 6, units = "in")

# very very top outliers
tab1_999 <- table(apply(aic_c_compare[zTz3 > quantile(zTz3, .999),], 1, which.min))/nrow(aic_c_compare[zTz3 > quantile(zTz3, .999),])
names(tab1_999) <- colnames(aic_c_compare)[as.numeric(unlist(dimnames(tab1_999)))]
png("plots/counts_model_win_999perc_outlier_loci_18models.png", width = 8, height = 6, res = 300, units = "in")
barplot(tab1_999, las = 2, main = "99.9% zTz top outlier loci - which model has lowest AICc?")
dev.off()

# than hmex or elevational trend (nothing too surprising here)
table(apply(aic_c_compare[zTz3 > quantile(zTz3, .99),], 1, which.min))/sum(zTz3 > quantile(zTz3, .99))
table(apply(aic_c_compare[zTz3 > quantile(zTz3, .95),], 1, which.min))/sum(zTz3 > quantile(zTz3, .95))
table(apply(aic_c_compare[zTz3 > quantile(zTz3, .90),], 1, which.min))/sum(zTz3 > quantile(zTz3, .9))
apply(aic_c_w_compare, 2, mean)
# for loci that seem to be selected, what are the weights across models? summarize across loci together to start
png("plots/mean_model_weights_99zTzoutlier.png", width = 8, height = 6, res = 300, units = "in")
model_weights <- apply(aic_c_w_compare[zTz3 > quantile(zTz3, .99),], 2, sum)/sum(zTz3 > quantile(zTz3, .99))
barplot(model_weights, las = 2, main = "mean model AICc weights across zTz>99% outlier loci")
dev.off()

# set 'other' category for loci with low bayes factors
m <- apply(aic_c_compare, 1, which.min)
table(m[sapply(1:length(m), function(i) bf[i, m[i]]) >= 100])/length(m)
sum(table(m[sapply(1:length(m), function(i) bf[i, m[i]]) >= 100])/length(m)) # about 13% of loci

table(m[sapply(1:length(m), function(i) bf[i, m[i]]) >= 1000])/length(m)
sum(table(m[sapply(1:length(m), function(i) bf[i, m[i]]) >= 1000])/length(m)) # about 3% of loci

# The above analyses compare relative fit but don't take into account MODEL ADEQUACY. 

# Can I show model predictions? I'll start by showing underlying raw data for loci that are top hits
bind_cols(maize_anc, sites) %>%
  bind_cols(., zb3_hmex) %>%
  mutate(., bf = bf$hmex) %>%
  filter(., apply(aic_c_compare, 1, which.min) == 3) %>%
  filter(., bf > 100) %>%
  .[c(T, rep(F, 100)), ] %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = mex_freq, shape = inv4m, color = zEnv, group = pos)) +
  geom_point() +
  geom_line() +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset loci where best model by AICc is hmex, bf > 100")

bind_cols(maize_anc, sites) %>%
  bind_cols(., zb3_elev) %>%
  mutate(., bf = bf$elev) %>%
  filter(., apply(aic_c_compare, 1, which.min) == 2) %>%
  filter(., bf > 100) %>%
  .[c(T, rep(F, 100)), ] %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = mex_freq, shape = inv4m, color = zEnv, group = pos)) +
  geom_point() +
  geom_line() +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset loci where best model by AICc is ~elev, bf > 100")

bind_cols(maize_anc, sites) %>%
  bind_cols(., zb3_elev) %>%
  filter(., apply(aic_c_compare, 1, which.min) == 2) %>%
  filter(., r_squared > .5) %>%
  .[c(T, rep(F, 10)), ] %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = mex_freq, shape = inv4m, color = zEnv, group = pos)) +
  geom_point() +
  geom_line() +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset loci where best model by AICc is ~elev, r^2 > .5")

# what do individual top outliers look like in unrotated space?

# one example locus from the inversion -- look at how well the model fits:
l1 <- data.frame(chr = 4, 
                 pos = 174081774)
l2 <- data.frame(chr = 4, pos = 45404997)
l3 <- data.frame(chr = 4, pos = 102315184)
for (l in c(l1, l2, l3)){
  a1 <- unlist(maize_anc[sites$chr == l["chr"] & sites$pos == l["pos"], ])

# fit elevation model
  zAnc1 = zAnc_maize$InvL %*% (a1 - zAnc_maize$alpha)
  zEnv1 = zAnc_maize$InvL %*% ordered_maize_pops$ELEVATION
  zInt1 = zAnc_maize$InvL %*% rep(1, 14)

  m1 <- lm(zAnc1 ~ 0 + zEnv1 + zInt1) # transformed intercept but no untransformed intercept
  fit1 <- data.frame(zAnc = zAnc1, 
                   zElev = zEnv1,
                   zInt = zInt1,
                   alpha = zAnc_maize$alpha,
                   zFit = m1$fitted.values) %>%
  bind_cols(., ordered_maize_pops)


#plot(fit1$alpha + fit1$zInt ~ fit1$ELEVATION)
#plot(t(chol(zAnc_maize$K))%*%fit1$zElev~fit1$ELEVATION) # transformation back works
  png(paste0("plots/model_pred_zElev_", l["chr"], "-", l["pos"], ".png"))
  plot(t(chol(zAnc_maize$K))%*%fit1$zFit+fit1$alpha~fit1$ELEVATION, 
     col = "blue", pch = 20, ylim = c(-.5, 1.5), 
     main = paste0("zAnc ~ zElev + zInt prediction vs. observed anc. zElev=", round(m1$coefficients[1],4), 
                   " zInt=", round(m1$coefficients[2],4)),
     ylab = "freq mex. ancestry", xlab = "elevation")
  points(a1~fit1$ELEVATION, col = "green", pch = 20)
  points(fit1$alpha~fit1$ELEVATION, col = "grey", pch = 20)
  abline(h = 0, col = "grey") # any predicted points outside of 0,1?
  abline(h = 1, col = "grey")
  legend(x = "topleft", 
       legend = c("pop mean ancestry", "observed ancestry", "model predicted ancestry"), 
       col = c("grey", "green", "blue"), 
       pch = 20)
  dev.off()
}
# total change in ancestry across our elevational gradient predicted by model:
max(t(chol(zAnc_maize$K))%*%fit1$zFit) - min(t(chol(zAnc_maize$K))%*%fit1$zFit) # 55% over ~ 1000 meters
# (55% beyond what can be explained by genomewide ancestry increasing with elevation)
# actually I don't want to include effects of the intercept..
max(a1) - min(a1) # measured change is ~85%

# TO DO: plot model expectations and true data for a large set of outliers
# how believable are my outliers and how big of a deal is it for estimates to > 1 or < 0
# I'll need to return fit+alpha for many model fits

# ~zElev individual loci
bind_cols(maize_anc, sites) %>%
  bind_cols(., zb3_elev) %>%
  filter(., pval_zEnv < .01) %>%
  .[c(T, rep(F, 100)), ] %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = mex_freq, shape = inv4m, color = zEnv, group = pos)) +
  geom_point() +
  geom_line() +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset low p-val outliers from zElev test - with rotated intercept zInt")
ggsave(paste0("plots/a_few_example_loci_zElev_outliers_rotated_intercept.png"),
       device = "png",
       width = 10, height = 8, units = "in")
# what is going on with high zEnv estimates but not very low pvalues?
# not great fits but still have steep slopes it looks like
bind_cols(maize_anc, sites) %>%
  bind_cols(., zb3_elev) %>%
  filter(., zEnv > .0005) %>%
  .[c(T, rep(F, 100)), ] %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = mex_freq, 
             shape = inv4m, color = log10(pval_zEnv), group = pos)) +
  geom_point() +
  geom_line() +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset largest pos. slopes zElev test - with rotated intercept zInt -show pvals")

# top pvals:
bind_cols(maize_anc, sites) %>%
  bind_cols(., zb3_elev) %>%
  filter(., pval_zEnv < .001) %>%
  .[c(T, rep(F, 10)), ] %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = mex_freq, 
             shape = inv4m, color = log10(pval_zEnv), group = pos)) +
  geom_point() +
  geom_line() +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset lowest p-vals zElev test - with rotated intercept zInt -show pvals")
# I think the issue is that steeper slopes aren't good fits because they're really logistic not linear
# but I don't know that's really true that they're logistic in transformed space.


# what do loci look like with higher 
bind_cols(maize_anc, sites) %>%
  bind_cols(., zb3_over1900m) %>%
  filter(., pval_zEnv < .001) %>%
  .[c(T, rep(F, 100)), ] %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = mex_freq, 
             shape = inv4m, color = zEnv, group = pos)) +
  geom_point() +
  geom_line() +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset low p-val outliers from over1900m test - with rotated intercept zInt")



# all selected
bind_cols(maize_anc, sites) %>%
  bind_cols(., zb3_hmex) %>%
  filter(., pval_zEnv < .01) %>%
  .[c(T, rep(F, 100)), ] %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = mex_freq, shape = inv4m, color = zEnv, group = pos)) +
  geom_point() +
  geom_line() +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset low p-val outliers from zAnc ~ all pops")
ggsave(paste0("plots/a_few_example_loci_zhmex_outliers.png"),
       device = "png",
       width = 10, height = 8, units = "in")

# one pop selected:
lapply(1:14, function(i){
  bind_cols(maize_anc, sites) %>%
    bind_cols(., zb3_onepop[[i]]) %>%
    filter(., r_squared > .5) %>% # extra filter to reduce noise
    filter(., pval_zEnv < .01) %>%
    .[c(T, rep(F, 100)), ] %>%
    gather(., "pop", "mex_freq", maize_pops) %>%
    left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
    ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = mex_freq, shape = inv4m, color = zEnv, group = pos)) +
    geom_point() +
    geom_line() +
    facet_wrap(~chr) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(paste("subset low p-val outliers from zAnc ~ one pop:", ordered_maize_pops$LOCALITY[i]))
  ggsave(paste0("plots/a_few_example_loci_onepop_outliers_", ordered_maize_pops$LOCALITY[i], ".png"),
         device = "png",
         width = 10, height = 8, units = "in")
}
  )







# plot outliers against their population means
d %>%
  bind_cols(., zb3_hmex) %>%
  filter(., pval_zEnv < .01) %>%
  ggplot(aes(x = pop_meanAlpha_maize, y = zEnv, color = log10(pval_zEnv))) +
  geom_point(size = .1) +
  ggtitle("zAnc outliers, p < .01")
ggsave(paste0("plots/zAnc_top_mex_outliers_in_maize_1percent.png"),
       device = "png",
       width = 10, height = 8, units = "in")
d %>%
  bind_cols(., zb3_hmex) %>%
  ggplot(aes(x = pop_meanAlpha_maize, y = zEnv, color = log10(pval_zEnv))) +
  geom_point(size = .1) +
  ggtitle("zAnc vs. population mean mexicana")
ggsave(paste0("plots/zAnc_pop_mean_mexicana_in_maize_allLoci.png"),
       device = "png",
       width = 10, height = 8, units = "in")

# outliers for zElev in and out of inversion
bind_cols(maize_anc, sites) %>%
  bind_cols(., zb3_elev) %>%
  filter(., pval_zEnv < .01) %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = mex_freq, fill = inv4m)) +
  geom_boxplot() +
  facet_wrap(~(zEnv >0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

d %>% # color by sum of squared residuals. okay actually that's only relevant as a relative value compared to the null model residuals
  bind_cols(., zb3_hmex) %>%
  filter(., r_squared > .5) %>%
  ggplot(aes(x = pop_meanAlpha_maize, y = zEnv, color = r_squared)) +
  geom_point(size = .1)



d %>%
  mutate(zTz = zTz3) %>%
  bind_cols(., zb3_hmex) %>%
  filter(zTz >= quantile(zTz, .95)) %>% # limit to top 1%
  mutate(log_pval = log10(pval_zEnv)) %>%
  ggplot(aes(x = meanAlpha_maize, y = zEnv, color = log_pval)) +
  geom_point(size = .1) +
  ggtitle("only showing top 5% zTz outliers; universal selection model")
d %>% # note: inv4m is not in top 1% of zTz but it is in top 5% of zTz
  mutate(zTz = zTz3) %>%
  bind_cols(., zb3_elev) %>%
  filter(zTz >= quantile(zTz, .95)) %>% # limit to top 1%
  mutate(log_pval = log10(pval_zEnv)) %>%
  ggplot(aes(x = meanAlpha_maize, y = zEnv, color = log_pval)) +
  #ggplot(aes(x = meanAlpha_maize, y = zEnv, color = inv4m)) +
  geom_point(size = .1) +
  ggtitle("only showing top 5% zTz outliers; elevation gradient model")


# what does zEnv look like genomewide for elevation?
d %>%
  mutate(zTz = zTz3) %>%
  bind_cols(., zb3_hmex) %>%
  ggplot(aes(x = pos, y = zEnv, color = zTz, shape = inv4m)) +
  geom_point(size = .1) +
  ggtitle("zElev outliers across the genome") +
  facet_wrap(~chr)
d %>%
  mutate(zTz = zTz3) %>%
  bind_cols(., zb3_hmex) %>%
  mutate(., log_pval = log(pval)) %>%
  ggplot(aes(x = pos, y = zEnv, color = zTz)) +
  geom_point(size = .1) +
  ggtitle("zAnc outliers across the genome") +
  facet_wrap(~chr)
ggsave("plots/zAnc_in_maize_pops_genomewide.png",
       height = 10, width = 14,
       units = "in", device = "png")
d %>%
  mutate(zTz = zTz3) %>%
  bind_cols(., zb3_elev) %>%
  mutate(., log_pval = log(pval)) %>%
  ggplot(aes(x = pos, y = zEnv, color = zTz)) +
  geom_point(size = .1) +
  ggtitle("zElev outliers across the genome") +
  facet_wrap(~chr)
ggsave("plots/zElev_in_maize_pops_genomewide.png",
       height = 10, width = 14,
       units = "in", device = "png")

# with intercept zInt for elevation:
d %>%
  mutate(zTz = zTz3) %>%
  bind_cols(., zb3_elev_int) %>%
  mutate(., log_pval = log(pval_zEnv)) %>%
  ggplot(aes(x = pos, y = zEnv, color = (inv4m | inv9e_or_d | inv9d_or_e))) +
  geom_point(size = .1) + 
  ggtitle("zElev outliers across the genome (w/ intercept)") +
  facet_wrap(~chr) # color inversion on chr9 too .. do I have the correct coordinates for these? is the length reasonable? are there others?
ggsave("plots/zElev_in_maize_pops_genomewide_w_intercept.png",
       height = 10, width = 14, 
       units = "in", device = "png")

# what is the qualitative diff of centering? zb2 vs. zb3:
plot(zb2$slope.zEnv ~ zb3_elev$slope.zEnv) # eh blob
cor(zb2$slope.zEnv, zb3_elev$slope.zEnv) # very weak correlation, or = .04
plot(zb2$pval_slope ~ zb3_elev$pval_slope) # black box (everywhere)
cor(zb2$pval_slope, zb3_elev$pval_slope) # slightly negative, cor -.06
plot(zb2$sum_sq_res ~ zb3_elev$sum_sq_res) # highly correlated
cor(zb2$sum_sq_res, zb3_elev$sum_sq_res) # cor = .96
cor(zb3_elev_int$zEnv, zb3_elev$zEnv) # poorly correlated
# plot zb3_elev_int$zEnv vs. simple beta elevation
plot(zb3_elev$zEnv[c(T, rep(F, 100))], simple_bElev_anc$envWeights[c(T, rep(F, 100))])
# slopes are similar but pvalues are not ..
plot(zb3_elev$pval_zEnv[c(T, rep(F, 100))], simple_bElev_anc$pval_Env[c(T, rep(F, 100))])
plot(zb3_elev$sum_sq_res[c(T, rep(F, 100))], simple_bElev_anc$sum_sq_res[c(T, rep(F, 100))])
# explanatory power is similar

# make MVN simulations for maize
n = 100000
seed = 500
set.seed(seed)
mvn_maize = mvrnorm(n=n, # create some MVN data
                     mu = zAnc_maize$alpha, 
                     Sigma = zAnc_maize$K,
                     empirical = F)
mvn_01_maize = mvn_maize # truncate at bounds [0, 1]
sum(mvn_01_maize<0)/(n*14) # about 5% less than 0
sum(mvn_01_maize>1)/(n*14) # and <.001% more than 1
mvn_01_maize[mvn_01_maize < 0] <- 0
mvn_01_maize[mvn_01_maize > 1] <- 1
hist(apply(mvn_maize, 1, mean)) # before truncation

mvn_01_maize_mean <- apply(mvn_01_maize, 1, mean)
hist(mvn_01_maize_mean) # after truncation
abline(v = quantile(mvn_01_maize_mean, c(.99, .999)), col = c("blue", "blue"))
summary(mvn_01_maize_mean)

maize_anc_mean <- apply(maize_anc, 1, mean)
summary(maize_anc_mean)
hist(maize_anc_mean)
abline(v = quantile(mvn_01_maize_mean, c(.99, .999)), col = c("blue", "darkblue"))
# slight enrichment but nothing crazy. What's my 5% false-discovery threshold?
# And how many loci pass it?
mean(maize_anc_mean > quantile(mvn_01_maize_mean, .99)) # 2x enriched
mean(maize_anc_mean > quantile(mvn_01_maize_mean, .999)) # 5-6x enriched
mean(maize_anc_mean > quantile(mvn_01_maize_mean, .9999)) # 20x enriched
fdr <- function(p, data, sims){
  obs <- sum(data > quantile(sims, p))
  null <- (1-p)*length(data)
  null/obs
}
sapply(c(.9, .99, .999, .9999, .99999), function(x) fdr(p = x, data = maize_anc_mean, sims = mvn_01_maize_mean))
fdr2 <- function(a, data, sims){# takes in an ancestry
  obs <- sum(data > a)
  null <- sum(sims > a)/length(sims)*length(data)
  null/obs
}
fdr2_low <- function(a, data, sims){# takes in an ancestry, test for low ancestry
  obs <- sum(data < a)
  null <- sum(sims < a)/length(sims)*length(data)
  null/obs
}

# set FDR thresholds for high mexicana ancestry in maize
if (rerun_all_models){
  test_anc <- seq(mean(maize_anc_mean), max(maize_anc_mean), length.out = 10000) # keep in range observed to not divide by 0
  test_fdr <- sapply(test_anc, function(x) fdr2(x, data = maize_anc_mean, sims = mvn_01_maize_mean))
  png("plots/power_curve_high_mex_in_maize_mvn_sim.png", 
    width = 8, height = 6, units = "in", res = 300)
  plot(test_anc, test_fdr, 
     main = "FDR high-mex outliers calculated based on MVN simulation of maize",
     pch = 20, cex = .1)
  abline(h=c(.1, .05,.01), col = c("skyblue", "orange", "red"))
  legend("topright", legend = c("FDR = 0.1", "FDR = 0.05", "FDR = 0.01"), col = c("skyblue", "orange", "red"), lty = 1)
  dev.off()

  # set FDR threshold for low mexicana ancestry in maize
  test_anc_low <- seq(min(maize_anc_mean), mean(maize_anc_mean), length.out = 10000) # keep in range observed to not divide by 0
  test_fdr_low <- sapply(test_anc_low, function(x) fdr2_low(x, data = maize_anc_mean, sims = mvn_01_maize_mean))
  png("plots/power_curve_low_mex_in_maize_mvn_sim.png", 
      width = 8, height = 6, units = "in", res = 300)
  plot(test_anc_low, test_fdr_low, # we have no power to detect low mexicana outliers individually 
       main = "FDR low-mex outliers calculated based on MVN simulation of maize",
       pch = 20, cex = .1)
  abline(h=c(.1, .05,.01), col = c("skyblue", "orange", "red"))
  legend("topright", legend = c("FDR = 0.1", "FDR = 0.05", "FDR = 0.01"), col = c("skyblue", "orange", "red"), lty = 1)
  dev.off()

# calculate FDR thresholds (approx)
  FDRs <- data.frame(thresholds = sapply(fdr_range, function(p) max(test_anc[test_fdr>p], na.rm = T)),
                     FDR = fdr_range, model = "high_mexicana_ancestry")      
  write.table(FDRs,
              "results/FDRs_high_mex.txt", quote = F, col.names = T, row.names = F, sep = "\t")
} else{
  FDRs <- read.table("results/FDRs_high_mex.txt", header = T, stringsAsFactors = F, sep = "\t")
}
# is seletion for mexicana ancestry a dominant mode of selection?
# what portion of zTz outliers meet 5% FDR threshold for high mex (vs some other mode of selection)



# plot maize ancestry across the genome colored by FDR
d_maize_mean_anc <- sites %>%
  mutate(mean_mexicana_anc_in_maize = maize_anc_mean)
d_maize_mean_anc$significance <- cut(d_maize_mean_anc$mean_mexicana_anc_in_maize,
                                     breaks = c(0, FDRs$thresholds, 1), include.lowest = F, 
                                     labels = c("n.s.", "FDR < 0.1", "FDR < 0.05", "FDR < 0.01"))
# add in mexicana anc in mexicana
#d_maize_mean_anc$mean_mexicana_anc_in_mexicana = apply(mexcana_anc, 2, mean)

ggplot(d_maize_mean_anc, aes(pos, mean_mexicana_anc_in_maize, color = significance)) +
  geom_point(size = .1) +
  scale_colour_manual(values = c("black", "skyblue", "orange", "red")) + 
  facet_wrap(~chr) +
  geom_abline(slope = 0, intercept = quantile(maize_anc_mean, .99), linetype = "dashed", color = "grey") +
  geom_abline(slope = 0, intercept = quantile(maize_anc_mean, .01), linetype = "dashed", color = "grey") +
  geom_abline(slope = 0, intercept = mean(maize_anc_mean), linetype = "dotted", color = "grey")
ggsave("plots/mean_mex_anc_in_maize_FDR_10_05_01_mvn_sims_whole_genome.png", device = "png",
       height = 8, width = 12, units = "in")

# how unusual are the ancestry slopes?
meta.pops.maize <- left_join(data.frame(pop = maize_pops, stringsAsFactors = F),
                             meta.pops, by = "pop")
maize_pops == meta.pops.maize$pop & maize_pops == colnames(maize_anc) # good
corr_elev_mvn_maize <- apply(mvn_01_maize, 1, function(x) cor(x, 
                                                            meta.pops.maize$ELEVATION))
# standard deviation is zero for some simulations where all pops have ancestry 0 -- so I set those to correlation = 0.
cor(rep(0,14), meta.pops.maize$ELEVATION)
corr_elev_mvn_maize[is.na(corr_elev_mvn_maize)] <- 0
corr_elev_maize <- apply(maize_anc, 1, function(x) cor(x, 
                                                          meta.pops.maize$ELEVATION)) 
hist(corr_elev_mvn_maize)

# set FDR thresholds for correlations with environment
test_elev_corr <- seq(0, max(corr_elev_maize), length.out = 1000) # keep in range observed to not divide by 0
test_fdr_elev_corr <- sapply(test_elev_corr, function(x) 
  fdr2(x, data = corr_elev_maize, sims = corr_elev_mvn_maize))
png("plots/power_curve_corrElev_in_maize_mvn_sim.png", 
    width = 8, height = 6, units = "in", res = 300)
plot(test_elev_corr, test_fdr_elev_corr, 
     main = "FDR correlation with elevation calculated based on MVN simulation of maize",
     pch = 20, cex = .1)
abline(h=c(.1, .05,.01), col = c("darkgreen", "orange", "red")) # I think what's happening is only the inversion
# comes out as really enriched .. and it's non-indep so I don't know if that's believable
# ok but overly conservative test because it's correcting for mean mexicana ancestry 
# that is already highly correlated (.84) with elevation
legend("topleft", legend = c("FDR = 0.1", "FDR = 0.05", "FDR = 0.01"), col = c("darkgreen", "orange", "red"), lty = 1)
dev.off()

# plot neg relationship between mexicana ancestry and environment
test_elev_corr_neg <- seq(min(corr_elev_maize), 0, length.out = 1000) # keep in range observed to not divide by 0
test_fdr_elev_corr_neg <- sapply(test_elev_corr_neg, function(x) 
  fdr2_low(x, data = corr_elev_maize, sims = corr_elev_mvn_maize))
png("plots/power_curve_neg_corrElev_in_maize_mvn_sim.png", 
    width = 8, height = 6, units = "in", res = 300)
plot(test_elev_corr_neg, test_fdr_elev_corr_neg, 
     main = "FDR neg. correlation with elevation calculated based on MVN simulation of maize",
     pch = 20, cex = .1)
abline(h=c(.1, .05,.01), col = c("darkgreen", "orange", "red")) # I think what's happening is only the inversion
# comes out as really enriched .. and it's non-indep so I don't know if that's believable
# ok but overly conservative test because it's correcting for mean mexicana ancestry 
# that is already highly correlated (.84) with elevation
legend("topleft", legend = c("FDR = 0.1", "FDR = 0.05", "FDR = 0.01"), col = c("darkgreen", "orange", "red"), lty = 1)
dev.off()

# I can still plot empirical 1% cutoffs for association with elevation
sites %>%
  mutate(corr_Elev = corr_elev_maize) %>%
  ggplot(., aes(pos, corr_Elev)) +
  geom_point(size = .1, color = "black") +
  facet_wrap(~chr) +
  geom_abline(slope = 0, intercept = quantile(corr_elev_maize, .99), linetype = "dashed", color = "grey") +
  geom_abline(slope = 0, intercept = quantile(corr_elev_maize, .01), linetype = "dashed", color = "grey") +
  geom_abline(slope = 0, intercept = mean(corr_elev_maize), linetype = "dotted", color = "grey")
ggsave("plots/mean_bElev_in_maize_FDR_empirical_1percent_whole_genome.png", device = "png",
       height = 8, width = 12, units = "in")



fdr(.99, data = corr_elev_maize, sims = corr_elev_mvn_maize)
quantile(corr_elev_mvn_maize, c(.99))


# actually I don't want the correlation, I want the slope of the linear model, Anc ~ elev, 
# because some loci have more variance in ancestry freq. than others

if (rerun_all_models){ # rerun or just reload results
  simple_bElev_anc <- apply(maize_anc, 
                            1, function(x)
                              simple_env_regression(ancFreq = x, envWeights = meta.pops.maize$ELEVATION)) %>%
    t(.) %>%
    as.data.frame(.)
  write.table(simple_bElev_anc,
              paste0("results/models/", PREFIX, "/maize/simple_bElev_anc.txt"),
              sep = "\t",
              quote = F,
              col.names = T, row.names = F)
  simple_bElev_mvn <- apply(mvn_01_maize, 
                            1, function(x)
                              simple_env_regression(ancFreq = x, envWeights = meta.pops.maize$ELEVATION)) %>%
    t(.) %>%
    as.data.frame(.) # check transposed!
  write.table(simple_bElev_mvn,
              paste0("results/models/", PREFIX, "/maize/simple_bElev_mvn_seed", seed, ".txt"),
              sep = "\t",
              quote = F,
              col.names = T, row.names = F)
}else{
  simple_bElev_anc <- read.table(paste0("results/models/", PREFIX, "/maize/simple_bElev_anc.txt"),
                                 sep = "\t",
                                 stringsAsFactors = F, header = T)
  simple_bElev_mvn <- read.table(paste0("results/models/", PREFIX, "/maize/simple_bElev_mvn_seed", seed, ".txt"),
              sep = "\t", stringsAsFactors = F, header = T)
}


summary(as.data.frame(simple_bElev_mvn)$envWeights)
summary(as.data.frame(simple_bElev_anc)$envWeights)

if (rerun_all_models){
  # get FDR thresholds
  test_simple_bElev <- seq(0, max(simple_bElev_anc$envWeights), length.out = 1000) # keep in range observed to not divide by 0
  test_fdr_simple_bElev <- sapply(test_simple_bElev, function(x) 
    fdr2(x, data = simple_bElev_anc$envWeights, sims = simple_bElev_mvn$envWeights))
  png("plots/power_curve_simple_bElev_in_maize_mvn_sim.png", 
      width = 8, height = 6, units = "in", res = 300)
  plot(test_simple_bElev, test_fdr_simple_bElev, 
       main = "FDR pos. slope simple Anc ~ elev calculated based on MVN simulation of maize",
       pch = 20, cex = .1)
  abline(h=c(.1, .05,.01), col = c("skyblue", "orange", "red")) # I think what's happening is only the inversion
  # comes out as really enriched .. and it's non-indep so I don't know if that's believable
  # ok but overly conservative test because it's correcting for mean mexicana ancestry 
  # that is already highly correlated (.84) with elevation
  legend("topright", legend = c("FDR = 0.1", "FDR = 0.05", "FDR = 0.01"), col = c("skyblue", "orange", "red"), lty = 1)
  dev.off()
  # neg slope with env
  test_simple_bElev_neg <- seq(min(simple_bElev_anc$envWeights), 0, length.out = 1000) # keep in range observed to not divide by 0
  test_fdr_simple_bElev_neg <- sapply(test_simple_bElev_neg, function(x) 
    fdr2_low(x, data = simple_bElev_anc$envWeights, sims = simple_bElev_mvn$envWeights))
  png("plots/power_curve_simple_bElev_neg_in_maize_mvn_sim.png", 
      width = 8, height = 6, units = "in", res = 300)
  plot(test_simple_bElev_neg, test_fdr_simple_bElev_neg, 
       main = "FDR neg. slope simple Anc ~ elev calculated based on MVN simulation of maize",
       pch = 20, cex = .1)
  abline(h=c(.1, .05,.01), col = c("skyblue", "orange", "red")) # I think what's happening is only the inversion
  # comes out as really enriched .. and it's non-indep so I don't know if that's believable
  # ok but overly conservative test because it's correcting for mean mexicana ancestry 
  # that is already highly correlated (.84) with elevation
  legend("topleft", legend = c("FDR = 0.1", "FDR = 0.05", "FDR = 0.01"), col = c("skyblue", "orange", "red"), lty = 1)
  dev.off()
  
  
  FDRs_simple_bElev_pos <- sapply(fdr_range, function(p) max(test_simple_bElev[test_fdr_simple_bElev>p], na.rm = T))
  FDRs_simple_bElev_neg <- sapply(fdr_range, function(p) max(test_simple_bElev_neg[test_fdr_simple_bElev_neg<p], na.rm = T))
  FDRs_simple_bElev <- data.frame(thresholds = c(FDRs_simple_bElev_neg, FDRs_simple_bElev_pos[3:1]), 
                                  FDR = c(rev(fdr_range), fdr_range), 
                                  sign = c(rep("neg", 3), rep("pos", 3))) %>%
    mutate(FDR_label = paste("FDR <", FDR))
  write.table(FDRs_simple_bElev,
              "results/FDRs_simple_bElev.txt", quote = F, col.names = T, row.names = F, sep = "\t")
} else{
  FDRs_simple_bElev <- read.table("results/FDRs_simple_bElev.txt", 
                                  header = T, stringsAsFactors = F, sep = "\t")
}

qqplot(simple_bElev_anc$envWeights, simple_bElev_mvn$envWeights)
abline(a = 0, b = 1, col = "blue")

# plot slopes with environment colored by FDR
d_simple_bElev <- sites %>%
  mutate(slope_elev = simple_bElev_anc$envWeights)
d_simple_bElev$significance <- cut(d_simple_bElev$slope_elev,
                                     breaks = c(-1, FDRs_simple_bElev$thresholds, 1), include.lowest = F, 
                                     labels = c("FDR < 0.01", "FDR < 0.05", "FDR < 0.1", "n.s.", "FDR < 0.1", "FDR < 0.05", "FDR < 0.01"))

ggplot(d_simple_bElev, 
       aes(pos, slope_elev, color = significance)) +
  geom_point(size = .1) +
  scale_colour_manual(values = c("red", "orange", "skyblue", "black")) + 
  facet_wrap(~chr, scales = "free_x") +
  geom_abline(slope = 0, intercept = quantile(d_simple_bElev$slope_elev, .99), linetype = "dashed", color = "grey") +
  #geom_abline(slope = 0, intercept = quantile(d_simple_bElev$slope_elev, .01), linetype = "dashed", color = "grey") +
  #geom_abline(slope = 0, intercept = mean(d_simple_bElev$slope_elev), linetype = "dotted", color = "grey")
ggsave("plots/simple_bElev_mex_anc_in_maize_FDR_10_05_01_mvn_sims_whole_genome.png", device = "png",
       height = 6, width = 10, units = "in")

# just chr9 QTL 
qtl9 <- read.table("../data/known_QTL/chr9_bin4_mhl1_locus_v4.bed",
                   header = F, sep = "\t")
colnames(qtl9) <- c("chr", "start", "end")
d_simple_bElev %>%
  filter(chr == 9) %>%
  filter(pos >= qtl9$start & pos <= qtl9$end) %>%
  ggplot(., 
       aes(pos, slope_elev, color = significance)) +
  geom_point(size = .1) +
  geom_abline(slope = 0, intercept = mean(d_simple_bElev$slope_elev), linetype = "dashed", color = "darkgrey") +
  scale_colour_manual(values = c("red", "orange", "skyblue", "darkgrey")) + 
  ylab("slope ancestry ~ elevation") +
  ggtitle("introgression at high elevation candidate near mhl1 - showing map bin 4")
ggsave("plots/simple_bElev_mex_anc_in_maize_chr9_mhl1_candidate_region.png",
       device = "png", height = 4, width = 6)

# just chr4 inversion
d_simple_bElev %>%
  filter(chr == 4) %>%
  filter(pos >= inv[inv$ID == "inv4m", "start"]-2.5*10^6 & pos <= inv[inv$ID == "inv4m", "end"]+1.5*10^6) %>%
  ggplot(., 
         aes(pos, slope_elev, color = significance)) +
  geom_point(size = .1) +
  geom_abline(slope = 0, intercept = mean(d_simple_bElev$slope_elev), linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = inv[inv$ID == "inv4m", "start"], color = "darkgrey") +
  geom_vline(xintercept = inv[inv$ID == "inv4m", "end"], color = "darkgrey") +
  scale_colour_manual(values = c("red", "orange", "skyblue", "darkgrey")) + 
  ylab("slope ancestry ~ elevation") +
  ggtitle("introgression across chr4 inversion")
ggsave("plots/simple_bElev_mex_anc_in_maize_inv4m.png",
       device = "png", height = 4, width = 6)


# manhattan style plot

# Compute chromosome size
sites_cum <- sites %>% 
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  # join to sites
  left_join(sites[, c("chr", "pos")], ., by = c("chr")) %>%
  mutate( pos_cum=pos+tot)

axis_spacing = sites_cum %>%
  group_by(chr) %>% 
  summarize(center=( max(pos_cum) + min(pos_cum) ) / 2,
            start = min(pos_cum),
            end = max(pos_cum))


# make wide format                       
d_simple_bElev %>%
  left_join(., sites_cum, by = c("chr", "pos")) %>%
              mutate(even_chr = (chr %% 2 == 0)) %>%
             ggplot(., aes(pos_cum, slope_elev, 
                           color = even_chr, show.legend = F)) +
  #geom_point_rast(size = .1) +
  geom_point(size = .1) +
  xlab("bp position on chromosome") +
  ylab("slope mexicana ancestry ~ elevation") +
  scale_colour_manual(values = c("grey", "darkgrey")) + 
  geom_abline(slope = 0, intercept = FDRs_simple_bElev$thresholds, linetype = "dashed", color = rep(c("red", "orange", "skyblue"), 2)) +
  geom_abline(slope = 0, intercept = mean(d_simple_bElev$slope_elev), color = "darkgrey", linetype = "dotted") +
  #geom_abline(data = FDRs_simple_bElev, aes(intercept = thresholds, slope = 0, color = FDR)) +
  #scale_colour_manual(values = c("grey", "darkgrey", "yellow", "red", "skyblue")) + 
  scale_x_continuous(label = axis_spacing$chr, breaks= axis_spacing$center) +
  theme(legend.position = "none")
ggsave("plots/simple_bElev_mex_anc_in_maize_FDR_10_05_01_mvn_sims_whole_genome_wide.png", device = "png",
       height = 3, width = 12, units = "in")


# hack to color outliers and by chromosome:
d_simple_bElev %>%
  left_join(., sites_cum, by = c("chr", "pos")) %>%
  mutate(color_by = ifelse(significance == "n.s." & (chr %% 2 == 0), 
                           "n.s. - even chrom", significance)) %>%
  ggplot(., aes(pos_cum, slope_elev, 
                color = color_by, show.legend = F)) +
  geom_point(size = .1) + 
  xlab("bp position on chromosome") +
  ylab("slope mexicana ancestry ~ elevation") +
  scale_colour_manual(values = c("red", "orange", "skyblue", "grey", "darkgrey")) + 
  #scale_colour_manual(values = c("grey", "darkgrey")) + 
  #geom_abline(slope = 0, intercept = FDRs_simple_bElev$thresholds, linetype = "dashed", alpha = .5,
  #            color = rep(c("red", "orange", "skyblue"), 2)) +
  geom_abline(slope = 0, intercept = mean(d_simple_bElev$slope_elev), color = "darkgrey", linetype = "dotted") +
  scale_x_continuous(label = axis_spacing$chr, breaks= axis_spacing$center) +
  theme(legend.position = "none")
ggsave("plots/simple_bElev_mex_anc_in_maize_FDR_10_05_01_mvn_sims_whole_genome_wide_color.png", device = "png",
       height = 3, width = 12, units = "in")
ggsave("plots/simple_bElev_mex_anc_in_maize_FDR_10_05_01_mvn_sims_whole_genome_less_wide_color.png", device = "png",
       height = 3, width = 12, units = "in")

# highlighting inv4m:
# where does chr4 start?
start_chr <- sites %>% 
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  group_by(chr) %>%
  summarise(start_chr = mean(tot)) %>%
  data.frame()

d_simple_bElev %>%
  left_join(., sites_cum, by = c("chr", "pos")) %>%
  mutate(color_by = ifelse(significance == "n.s." & (chr %% 2 == 0), 
                           "n.s. - even chrom", significance)) %>%
  ggplot(., aes(pos_cum, slope_elev, 
                color = color_by, show.legend = F)) +
  geom_point(size = .1) +
  xlab("bp position on chromosome") +
  ylab("slope mexicana ancestry ~ elevation") +
  scale_colour_manual(values = c("red", "orange", "skyblue", "grey", "darkgrey")) + 
  geom_vline(xintercept = inv[inv$ID == "inv4m", "start"] + start_chr[start_chr$chr==4, "start_chr"], color = "black", linetype = "dashed") +
  geom_vline(xintercept = inv[inv$ID == "inv4m", "end"] + start_chr[start_chr$chr==4, "start_chr"], color = "black", linetype = "dashed") +
  geom_abline(slope = 0, intercept = mean(d_simple_bElev$slope_elev), color = "darkgrey", linetype = "dotted") +
  scale_x_continuous(label = axis_spacing$chr, breaks= axis_spacing$center) +
  theme(legend.position = "none")
ggsave("plots/simple_bElev_mex_anc_in_maize_FDR_10_05_01_mvn_sims_inv4m_wide_color.png", device = "png",
       height = 3, width = 12, units = "in")
ggsave("plots/simple_bElev_mex_anc_in_maize_FDR_10_05_01_mvn_sims_inv4m_less_wide_color.png", device = "png",
       height = 3, width = 9, units = "in")
# highlighting chr9 mhl1 locus
d_simple_bElev %>%
  left_join(., sites_cum, by = c("chr", "pos")) %>%
  mutate(color_by = ifelse(significance == "n.s." & (chr %% 2 == 0), 
                           "n.s. - even chrom", significance)) %>%
  ggplot(., aes(pos_cum, slope_elev, 
                color = color_by, show.legend = F)) +
  geom_point(size = .1) +
  xlab("bp position on chromosome") +
  ylab("slope mexicana ancestry ~ elevation") +
  scale_colour_manual(values = c("red", "orange", "skyblue", "grey", "darkgrey")) + 
  geom_vline(xintercept = qtl9$start + start_chr[start_chr$chr==9, "start_chr"], color = "black", linetype = "dashed") +
  geom_vline(xintercept = qtl9$end + start_chr[start_chr$chr==9, "start_chr"], color = "black", linetype = "dashed") +
  geom_abline(slope = 0, intercept = mean(d_simple_bElev$slope_elev), color = "darkgrey", linetype = "dotted") +
  scale_x_continuous(label = axis_spacing$chr, breaks= axis_spacing$center) +
  theme(legend.position = "none")
ggsave("plots/simple_bElev_mex_anc_in_maize_FDR_10_05_01_mvn_sims_mhl1_wide_color.png", device = "png",
       height = 3, width = 12, units = "in")
ggsave("plots/simple_bElev_mex_anc_in_maize_FDR_10_05_01_mvn_sims_mhl1_less_wide_color.png", device = "png",
       height = 3, width = 9, units = "in")

  
  
  
d_maize_mean_anc %>% # note no sig on low side because simulations exceeded values found in real data (namely, 0)
  left_join(., sites_cum, by = c("chr", "pos")) %>%
  mutate(even_chr = (chr %% 2 == 0)) %>%
  ggplot(., aes(pos_cum, mean_mexicana_anc_in_maize, 
                color = even_chr, show.legend = F)) +
  geom_point(size = .1) +
  xlab("bp position on chromosome") +
  ylab("mexicana ancestry frequency") +
  ylim(c(0,1)) +
  scale_colour_manual(values = c("grey", "darkgrey")) + 
  geom_abline(slope = 0, intercept = c(FDRs), linetype = "dashed", color = c("skyblue", "orange", "red")) +
  geom_abline(slope = 0, intercept = mean(d_maize_mean_anc$mean_mexicana_anc_in_maize), color = "darkgrey", linetype = "dotted") +
  scale_x_continuous(label = axis_spacing$chr, breaks= axis_spacing$center) +
  theme(legend.position = "none")
ggsave("plots/mean_mex_anc_in_maize_FDR_10_05_01_mvn_sims_whole_genome_wide.png", device = "png",
       height = 3, width = 12, units = "in")
mean(simple_bElev_anc$envWeights)
# hack to color outliers
d_maize_mean_anc %>% # note no sig on low side because simulations exceeded values found in real data (namely, 0-1)
  left_join(., sites_cum, by = c("chr", "pos")) %>%
  mutate(color_by = ifelse(significance == "n.s." & (chr %% 2 == 0), 
                           "n.s. - even chrom", significance)) %>%
  ggplot(., aes(pos_cum, mean_mexicana_anc_in_maize, 
                color = color_by, show.legend = F)) +
  geom_point(size = .1) +
  xlab("bp position on chromosome") +
  ylab("mexicana ancestry frequency") +
  ylim(c(0,1)) +
  scale_colour_manual(values = c("grey", "skyblue", "orange", "red", "darkgrey")) + 
  #scale_colour_manual(values = c("grey", "darkgrey")) + 
  #geom_abline(slope = 0, intercept = c(FDRs), linetype = "dashed", alpha = .5, 
  #            color = c("skyblue", "orange", "red")) +
  geom_abline(slope = 0, intercept = mean(d_maize_mean_anc$mean_mexicana_anc_in_maize), color = "darkgrey", linetype = "dotted") +
  scale_x_continuous(label = axis_spacing$chr, breaks= axis_spacing$center) +
  theme(legend.position = "none")
ggsave("plots/mean_mex_anc_in_maize_FDR_10_05_01_mvn_sims_whole_genome_wide_color.png", device = "png",
      height = 3, width = 12, units = "in")
ggsave("plots/mean_mex_anc_in_maize_FDR_10_05_01_mvn_sims_whole_genome_less_wide_color.png", device = "png",
       height = 3, width = 9, units = "in")


# look at pvalues for bElev instead of bElev directly:
test_simple_bElev_log10_pval <- seq(0, max(-log10(simple_bElev_anc$pval_Env)), length.out = 1000) # keep in range observed to not divide by 0
test_fdr_simple_bElev_log10_pval <- sapply(test_simple_bElev_log10_pval, function(x) 
  fdr2(x, data = -log10(simple_bElev_anc$pval_Env), 
           sims = -log10(simple_bElev_mvn$pval_Env)))
#png("plots/power_curve_simple_bElev_pval_in_maize_mvn_sim.png", 
#    width = 8, height = 6, units = "in", res = 300)
plot(test_simple_bElev_log10_pval, test_fdr_simple_bElev_log10_pval, 
     main = "FDR log10 pval Anc ~ elev calculated based on MVN simulation of maize",
     pch = 20, cex = .1)
abline(h=c(.1, .05,.01), col = c("skyblue", "orange", "red")) # I think what's happening is only the inversion
#dev.off()
png("plots/pval_vs_slope_simple_bElev_in_maize_data_and_mvn_sim.png", 
        width = 8, height = 6, units = "in", res = 300)
plot(simple_bElev_anc$envWeights, -log10(simple_bElev_anc$pval_Env), 
     col = alpha("orange", .1),
     main = "OLS pval vs. slope mexicana anc. with elevation",
     sub = "orange = data, grey = simulations",
     xlab = "OLS slope w/ elevation",
     ylab = "OLS pval for slope")
points(simple_bElev_mvn$envWeights, -log10(simple_bElev_mvn$pval_Env), col = alpha("grey", .1))
dev.off() # ok to use sims as null. these pvals aren't very meaningful given the independence assumption isn't met.
# better to go with the simulated data than the parametric approximation pval

# make sure the simulated MVN has the expected chi-squared distribution for zTz ~ chi-squared(14)
# chisq expectation
chisq14 <- rchisq(n=100000, df=14)

# calculate zAnc for mvn simulation
zAnc3_mvn <- t(apply(mvn_maize, 1, function(l) zAnc(ancFreq = l, 
                                                          invL = zAnc_maize$InvL, 
                                                          alpha = zAnc_maize$alpha)))
# zTz
zTz3_mvn <- apply(zAnc3_mvn, 1, function(i) t(i) %*% i)
hist(zTz3_mvn)
hist(chisq14, add = T, col = "blue") # fits untruncated mvn

# truncated mvn
zAnc3_01_mvn <- t(apply(mvn_01_maize, 1, function(l) zAnc(ancFreq = l, 
                                                          invL = zAnc_maize$InvL, 
                                                          alpha = zAnc_maize$alpha)))

zTz3_01_mvn <- apply(zAnc3_01_mvn, 1, function(i) t(i) %*% i)
hist(zTz3_01_mvn, freq = T)
hist(chisq14, add = T, col = alpha("blue", .3)) # fits untruncated mvn



quantile(chisq14, c(.9, .99, .999, .9999))
quantile(zTz3_mvn, c(.9, .99, .999, .9999))
quantile(zTz3_01_mvn, c(.9, .99, .999, .9999))
sum(zTz3_01_mvn > quantile(chisq14, .99))/length(zTz3_01_mvn) # truncation reduces outliers slightly
# so the result is a slightly conservative test.

# what does my data look like?
# compared to these others?
png("plots/zTz_chisq14_hist_dist_MVN_sim_and_data.png", 
    width = 4, height = 12, units = "in", res = 300)
par(mfrow = c(3, 1))
hist(zTz3_mvn, freq = F, ylim = c(0, .1), breaks = 30, main = "zTz MVN sim. - no truncation")
curve(expr = dchisq(x, df = 14), from = 0, to = 100,
      col = "blue", add = T)

hist(zTz3_01_mvn, freq = F, ylim = c(0, .1), breaks = 30, "zTz MVN sim. - [0, 1] truncation")
curve(expr = dchisq(x, df = 14), from = 0, to = 100,
      col = "blue", add = T)

hist(zTz3, freq = F, ylim = c(0, .1), breaks= 30, "zTz maize data")
curve(expr = dchisq(x, df = 14), from = 0, to = 100,
      col = "blue", add = T)
par(mfrow = c(1, 1))
dev.off()

# ok so using zTz outliers, how many of them appear to fit
# a steep elevational gradient?
sum(zTz3 > quantile(chisq14, .999))/length(zTz3)
# fdr at .999 threshold? about 2%
.001*length(zTz3)/(sum(zTz3 > quantile(chisq14, .999)))

test_zTz <- seq(0, max(zTz3), length.out = 1000) # keep in range observed to not divide by 0
test_fdr_zTz <- sapply(test_zTz, 
                       function(x) (1 - pchisq(q = x, df = 14))*length(zTz3)/(sum(zTz3 > x))) 
png("plots/power_curve_zTz_chisq14.png", 
    width = 8, height = 6, units = "in", res = 300)
plot(test_zTz, test_fdr_zTz, 
     main = "FDR zTz outliers calculated based on chi-squared dist. df=14",
     pch = 20, cex = .1)
abline(h=c(.1, .05,.01), col = c("skyblue", "orange", "red"))
dev.off()

# calculate FDR thresholds (approx)
FDRs_zTz_chisq14_approx <- data.frame(thresholds = sapply(fdr_range, function(p) min(test_zTz[test_fdr_zTz<p], na.rm = T)),
                                      FDR = fdr_range, model = "zTz_chisq14")
write.table(FDRs_zTz_chisq14_approx,
            "results/FDRs_zTz_chisq14.txt", quote = F, col.names = T, row.names = F, sep = "\t")
FDRs_zTz_chisq14_approx <- read.table("results/FDRs_zTz_chisq14.txt", header = T, sep = "\t", stringsAsFactors = F)


# plot zTz across the maize genome:
data.frame(sites, zTz = zTz3) %>%
  left_join(., sites_cum, by = c("chr", "pos")) %>%
  mutate(even_chr = (chr %% 2 == 0)) %>%
  ggplot(., aes(pos_cum, zTz, 
                color = even_chr, show.legend = F)) +
  geom_point(size = .1) +
  xlab("bp position on chromosome") +
  ylab("zTz statistic") +
  scale_colour_manual(values = c("grey", "darkgrey")) + 
  geom_abline(slope = 0, intercept = FDRs_zTz_chisq14_approx, linetype = "dashed", color = c("red", "orange", "skyblue")) +
  scale_x_continuous(label = axis_spacing$chr, breaks= axis_spacing$center) +
  theme(legend.position = "none")
ggsave("plots/zTz_in_maize_FDR_10_05_01_using_chisq14_whole_genome_wide.png", device = "png",
       height = 3, width = 12, units = "in")

# draw qqplot
png("plots/qqplot_zTz_mvn_01_chisq14.png", 
    width = 8, height = 6, units = "in", res = 300)
qqplot(x = qchisq(ppoints(500), df = 14), y = zTz3_01_mvn,
       main = "QQ plot for truncated [0,1] mvn simulation zTz statistic vs. chi-sq df =14")
qqline(y = zTz3_01_mvn, distribution = function(p) qchisq(p, df = 14),
       prob = c(0.25, 0.75), col = "blue")
#abline(a = 0, b = 1, col = "blue")
dev.off()
png("plots/qqplot_zTz_mvn_01_vs_maize_data.png", 
    width = 8, height = 6, units = "in", res = 300)
qqplot(x = zTz3_01_mvn, y = zTz3, main = "QQ plot for observed zTz statistic vs. simulated mvn expectation",
       sub = "line = chi-sq df=14 expectation")
qqline(y = zTz3, distribution = function(p) qchisq(p, df = 14),
       prob = c(0.25, 0.75), col = "blue")
abline(a = 0, b = 1, col = "blue")
dev.off()
table(zTz3 > 20)/length(zTz3)
table(zTz3 > 30)/length(zTz3)

# can I say what % of the top 1% FDR in zTz overlap with top 1% fdr for high mexicana ancestry or for anc ~ elev?
top1_outliers <- sum(zTz3 > FDRs_zTz_chisq14_approx[1])
table(zTz3 > filter(FDRs_zTz_chisq14_approx, FDR == 0.01)[["thresholds"]], maize_anc_mean > filter(FDRs, FDR == 0.01)[["thresholds"]])
table(zTz3 > filter(FDRs_zTz_chisq14_approx, FDR == 0.01)[["thresholds"]], simple_bElev_anc$envWeights < filter(FDRs_simple_bElev, sign == "neg" & FDR == 0.01)[["thresholds"]])
table(zTz3 > filter(FDRs_zTz_chisq14_approx, FDR == 0.01)[["thresholds"]], simple_bElev_anc$envWeights > filter(FDRs_simple_bElev, sign == "pos" & FDR == 0.01)[["thresholds"]])
table(zTz3 > filter(FDRs_zTz_chisq14_approx, FDR == 0.01)[["thresholds"]], simple_bElev_anc$envWeights < filter(FDRs_simple_bElev, sign == "neg" & FDR == 0.01)[["thresholds"]], maize_anc_mean > FDRs[1])

# top 1% FDR outliers
outliers_top01 <- data.frame(zTz = zTz3 > filter(FDRs_zTz_chisq14_approx, FDR == 0.01)[["thresholds"]],
                       pos_Elev = simple_bElev_anc$envWeights > filter(FDRs_simple_bElev, sign == "pos" & FDR == 0.01)[["thresholds"]],
                       neg_Elev = simple_bElev_anc$envWeights < filter(FDRs_simple_bElev, sign == "neg" & FDR == 0.01)[["thresholds"]],
                       high_mex = maize_anc_mean > filter(FDRs, FDR == 0.01)[["thresholds"]]) %>%
  mutate(both = high_mex & pos_Elev) %>%
  filter(zTz) %>%
  apply(., 2, sum) %>%
  t(.) %>%
  .[1,]
outliers_top01 %>%
  write.table(., "results/outlier_overlap_zTz.txt", quote = F, col.names = T, row.names = F, sep = "\t")

outliers_top05 <- data.frame(zTz = zTz3 > filter(FDRs_zTz_chisq14_approx, FDR == 0.05)[["thresholds"]],
                             pos_Elev = simple_bElev_anc$envWeights > filter(FDRs_simple_bElev, sign == "pos" & FDR == 0.05)[["thresholds"]],
                             neg_Elev = simple_bElev_anc$envWeights < filter(FDRs_simple_bElev, sign == "neg" & FDR == 0.05)[["thresholds"]],
                             high_mex = maize_anc_mean > filter(FDRs, FDR == 0.05)[["thresholds"]]) %>%
  filter(zTz) %>%
  mutate(both = high_mex & pos_Elev) %>%
  apply(., 2, sum) %>%
  t(.) %>%
  .[1,]
outliers_top05 %>%
  write.table(., "results/outlier_overlap_zTz_FRD5.txt", quote = F, col.names = T, row.names = F, sep = "\t")



# like phylogenetic least squares

# plot some of the top hits:
bind_cols(d_simple_bElev) %>%
  bind_cols(maize_anc) %>%
  .[c(T, rep(F, 500)), ] %>%
  filter(inv4m) %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = ELEVATION, y = mex_freq, shape = inv4m, 
             color = significance, group = pos)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("red", "orange", "skyblue", "darkgrey")) +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset outliers slope anc ~ elevation")
ggsave(paste0("plots/a_few_example_loci_bElev_outliers_inv4m.png"),
       device = "png",
       width = 6, height = 4, units = "in")
# more outlier:
# plot some of the top hits:
bind_cols(d_simple_bElev) %>%
  bind_cols(maize_anc) %>%
  .[c(T, rep(F, 100)), ] %>%
  filter(inv4m) %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = ELEVATION, y = mex_freq, shape = inv4m, 
             color = significance, group = pos)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("red", "orange", "skyblue", "darkgrey")) +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset outliers slope anc ~ elevation")
ggsave(paste0("plots/a_few_more_example_loci_bElev_outliers_inv4m.png"),
       device = "png",
       width = 6, height = 4, units = "in")
# max vs.  genomewide background
bind_cols(d_simple_bElev) %>%
  bind_cols(maize_anc) %>%
  filter(slope_elev == max(slope_elev)) %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = ELEVATION, y = mex_freq, shape = inv4m, 
             group = pos, color = significance)) +
  geom_point() +
  geom_line() +
  geom_point(aes(y = zAnc_maize$alpha, x = arrange(meta.pops.maize, popN)$ELEVATION), color = "darkgrey") +
  geom_line(aes(y = zAnc_maize$alpha, x = arrange(meta.pops.maize, popN)$ELEVATION), linetype = "dashed", color = "darkgrey") +
  scale_color_manual(values = c("red", "orange", "skyblue", "darkgrey")) +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset outliers slope anc ~ elevation")
ggsave(paste0("plots/top_outlier_bElev_outliers_inv4m_vs_genomewide.png"),
       device = "png",
       width = 6, height = 4, units = "in")


# plot ancestry in mhl1 region
bind_cols(d_simple_bElev) %>%
  bind_cols(maize_anc) %>%
  .[c(T, rep(F, 100)), ] %>%
  filter(chr==9 & pos > 109.5*10^6 & pos < 110.5*10^6) %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = mex_freq, shape = inv4m, 
             color = significance, group = pos)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("red", "orange", "skyblue", "darkgrey")) +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset outliers slope anc ~ elevation")
ggsave(paste0("plots/a_few_more_example_loci_bElev_outliers_mhl1_pops.png"),
       device = "png",
       width = 6, height = 4, units = "in")

bind_cols(d_simple_bElev) %>%
  bind_cols(maize_anc) %>%
  .[c(T, rep(F, 100)), ] %>%
  filter(chr==9 & pos > 109.5*10^6 & pos < 110.5*10^6) %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = ELEVATION, y = mex_freq, shape = inv4m, 
             color = significance, group = pos)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("red", "orange", "skyblue", "darkgrey")) +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset outliers slope anc ~ elevation")
ggsave(paste0("plots/a_few_more_example_loci_bElev_outliers_mhl1.png"),
       device = "png",
       width = 6, height = 4, units = "in")

d %>%
  filter(meanAlpha_maize == max(meanAlpha_maize)) %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = ELEVATION, y = mex_freq, shape = inv4m, 
             group = pos)) +
  geom_point(color = "red") +
  geom_line(color = "red") +
  geom_point(aes(y = zAnc_maize$alpha, x = arrange(meta.pops.maize, popN)$ELEVATION), color = "darkgrey") +
  geom_line(aes(y = zAnc_maize$alpha, x = arrange(meta.pops.maize, popN)$ELEVATION), linetype = "dashed", color = "darkgrey") +
  scale_color_manual(values = c("red", "orange", "skyblue", "darkgrey")) +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset outliers slope anc ~ elevation")
ggsave(paste0("plots/top_outlier_loci_high_mexicana.png"),
       device = "png",
       width = 6, height = 4, units = "in")

d %>%
  filter(meanAlpha_maize == max(meanAlpha_maize)) %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = ELEVATION, y = mex_freq, shape = inv4m, 
             group = pos)) +
  geom_point(color = "red") +
  geom_line(color = "red") +
  geom_point(aes(y = zAnc_maize$alpha, x = arrange(meta.pops.maize, popN)$ELEVATION), color = "darkgrey") +
  geom_line(aes(y = zAnc_maize$alpha, x = arrange(meta.pops.maize, popN)$ELEVATION), linetype = "dashed", color = "darkgrey") +
  scale_color_manual(values = c("red", "orange", "skyblue", "darkgrey")) +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset outliers slope anc ~ elevation")
ggsave(paste0("plots/top_outlier_loci_high_mexicana.png"),
       device = "png",
       width = 6, height = 4, units = "in")
# plotting the numbers I got out of the model for top hits:
top_ztz_hits <- data.frame(hit = 1:outliers_top01[["zTz"]],
                           type = c(rep("high mexicana", outliers_top01[["high_mex"]]),
                                    rep("higher mexicana w/ elev.", outliers_top01[["pos_Elev"]]),
                                    rep("lower mexicana w/ elev.", outliers_top01[["neg_Elev"]]),
                                    rep("Other", (outliers_top01[["zTz"]] - sum(outliers_top01[c("high_mex", "pos_Elev", "neg_Elev")])))))
ggplot(top_ztz_hits, aes(x = type, fill = type)) +
  #geom_bar(stat = "percent") +
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("relative frequencies in zTz outliers") +
  ggtitle("overlap of top 1% FDR zTz hits and other outliers") +
  theme(axis.text.x = element_text(angle = 90))
ggsave("plots/overlap_barplot_ztz_top_and_other_sel_models_top.png",
       device = "png", height = 6, width = 6)
outliers_top01[["zTz"]]/length(zTz3) #2.7%

top_ztz_hits05 <- data.frame(hit = 1:outliers_top05[["zTz"]],
                           type = c(rep("high mexicana", outliers_top05[["high_mex"]]),
                                    rep("higher mexicana w/ elev.", outliers_top05[["pos_Elev"]]),
                                    rep("lower mexicana w/ elev.", outliers_top05[["neg_Elev"]]),
                                    rep("Other", (outliers_top05[["zTz"]] - sum(outliers_top05[c("high_mex", "pos_Elev", "neg_Elev")])))))
ggplot(top_ztz_hits05, aes(x = type, fill = type)) +
  #geom_bar(stat = "percent") +
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("relative frequencies in zTz outliers") +
  ggtitle("overlap of top 5% FDR zTz hits and other outliers") +
  theme(axis.text.x = element_text(angle = 90))
ggsave("plots/overlap_barplot_ztz_top_and_other_sel_models_top_FDR5.png",
       device = "png", height = 6, width = 6)




flwr_hits <- data.frame(
  chr6_centromere = 52*10^6 + axis_spacing$start[axis_spacing$chr == 6],
chr3_centromere = 86*10^6 + axis_spacing$start[axis_spacing$chr == 3],
chr3_inversion = 79*10^6 + axis_spacing$start[axis_spacing$chr == 3],
chr5_centromere = 105*10^6 + axis_spacing$start[axis_spacing$chr == 5],
chr4_inv_start = 171771502 + axis_spacing$start[axis_spacing$chr == 4],
chr4_inv_end = 185951149 + axis_spacing$start[axis_spacing$chr == 4]) %>%
  tidyr::gather(., "sv", "pos_cum")
chrom_positions <- axis_spacing %>%
  dplyr::select(c("chr", "start", "end"))%>%
  tidyr::gather(., "name", "pos_cum", 2:3) %>%
  mutate(sv = paste0(chr, "_", name)) %>%
  bind_rows(., flwr_hits)

chrom_positions %>%
  ggplot(., aes(x = pos_cum, y = 1, color = sv)) +
  geom_point()
ggsave("plots/flower_time_structural_variants.png",
       device = "png", height = 3, width = 12)


# chr6 centromere 52Mb; 29%
# chr3 86Mb
# chr5 105Mb
# chr3 v3. 79Mb - 6Mb inversion??

# QQ plots
png("plots/QQ_mean_maize_anc_MVN01_data.png",
    height = 6, width = 6, units = "in", res = 300)
qqplot(mvn_01_maize_mean, maize_anc_mean, cex = .1,
       main = "QQ plot mean mexicana ancestry across maize pops",
       xlab = "mean mex. ancestry in maize MVN (truncated) sims",
       ylab = "mean mex. ancestry in maize from data")
abline(a = 0, b = 1, col = "blue")
#quantile(mvn_01_maize_mean, c(.5, .75, .8, .9, .95, .99, .999))
#quantile(maize_anc_mean, c(.5, .75, .8, .9, .95, .99, .999))
abline(v = quantile(mvn_01_maize_mean, c(.01, .05, .1, .9, .95, .99)),
       col = c("red", "orange", "skyblue", "skyblue", "orange", "red"))
legend("topleft", legend = c(.9, .95, .99), col = c("skyblue", "orange", "red"),
       lty = 1, title = "percentile")
dev.off()
# difference of truncating or not
png("plots/QQ_mean_maize_anc_MVN01_MVN.png",
    height = 6, width = 6, units = "in", res = 300)
qqplot(mvn_01_maize_mean, apply(mvn_maize, 1, mean), cex = .1,
       main = "QQ plot mean mexicana ancestry across maize pops")
abline(a = 0, b = 1, col = "blue")
#quantile(mvn_01_maize_mean, c(.5, .75, .8, .9, .95, .99, .999))
#quantile(maize_anc_mean, c(.5, .75, .8, .9, .95, .99, .999))
abline(v = quantile(apply(mvn_maize, 1, mean), c(.01, .05, .1, .9, .95, .99)),
       col = c("red", "orange", "skyblue", "skyblue", "orange", "red"))
legend("topleft", legend = c(.9, .95, .99), col = c("skyblue", "orange", "red"),
       lty = 1, title = "percentile")
dev.off()

# slope simple bElev:
png("plots/QQ_slope_bElev_MVN01_data.png",
    height = 6, width = 6, units = "in", res = 300)
qqplot(simple_bElev_mvn$envWeights, simple_bElev_anc$envWeights, cex = .1,
       main = "QQ plot slope mex anc ~ elev across maize pops",
       xlab = "slopes (truncated) MVN simulation",
       ylab = "slopes observed in data")
abline(a = 0, b = 1, col = "blue")
abline(v = quantile(simple_bElev_mvn$envWeights, c(.01, .05, .1, .9, .95, .99)),
       col = c("red", "orange", "skyblue", "skyblue", "orange", "red"))
legend("topleft", legend = c(.9, .95, .99), col = c("skyblue", "orange", "red"),
       lty = 1, title = "percentile")
dev.off()

log10_pvals_simple_bElev_mvn <- log10(simple_bElev_mvn$pval_Env)
log10_pvals_simple_bElev_mvn[is.na(log10_pvals_simple_bElev_mvn)] <- 0 # some are NA for slope zero

png("plots/QQ_pval_bElev_MVN01_data.png",
    height = 6, width = 6, units = "in", res = 300)
qqplot(log10_pvals_simple_bElev_mvn, log10(simple_bElev_anc$pval_Env), cex = .1,
       main = "QQ plot pval mex anc ~ elev across maize pops",
       xlab = "log10 pval slope w/ elev (truncated) MVN sim",
       ylab = "log10 pval slope w/ elev from data")
abline(a = 0, b = 1, col = "blue")
abline(v = quantile(log10_pvals_simple_bElev_mvn, c(.01, .05, .1, .9, .95, .99)),
       col = c("red", "orange", "skyblue", "skyblue", "orange", "red"))
legend("topleft", legend = c(.9, .95, .99), col = c("skyblue", "orange", "red"),
       lty = 1, title = "percentile")
dev.off()

# zTz QQ plots
png("plots/QQ_zTz_MVN01_chisq.png",
    height = 6, width = 6, units = "in", res = 300)
qqplot(#x = qchisq(ppoints(500), df = 14), 
      x = rchisq(1000000, df = 14),
       y = zTz3_01_mvn,
       main = "QQ plot truncated [0,1] mvn sim. zTz vs. chi-sq df =14",
       cex = .1)
#qqline(y = zTz3_01_mvn, distribution = function(p) qchisq(p, df = 14),
#       prob = c(0.25, 0.75), col = "blue")
abline(a = 0, b = 1, col = "blue")
abline(v = quantile(rchisq(100000, df = 14), c(.01, .05, .1, .9, .95, .99)),
       col = c("red", "orange", "skyblue", "skyblue", "orange", "red"))
legend("topleft", legend = c(.9, .95, .99), col = c("skyblue", "orange", "red"),
       lty = 1, title = "percentile")
dev.off()
# non-truncated MVN simulation vs. chisq:
png("plots/QQ_zTz_MVN_chisq.png",
    height = 6, width = 6, units = "in", res = 300)
qqplot(#x = qchisq(ppoints(500), df = 14), 
  x = rchisq(1000000, df = 14),
  y = zTz3_mvn,
  main = "QQ plot untruncated mvn sim. zTz vs. chi-sq df =14",
  cex = .1)
abline(a = 0, b = 1, col = "blue")
abline(v = quantile(rchisq(100000, df = 14), c(.01, .05, .1, .9, .95, .99)),
       col = c("red", "orange", "skyblue", "skyblue", "orange", "red"))
legend("topleft", legend = c(.9, .95, .99), col = c("skyblue", "orange", "red"),
       lty = 1, title = "percentile")
dev.off()

# QQ plot zTz data vs. chi-sq
png("plots/QQ_zTz_data_chisq.png",
    height = 6, width = 6, units = "in", res = 300)
qqplot(#x = qchisq(ppoints(500), df = 14), 
  x = rchisq(1000000, df = 14),
  y = zTz3,
  main = "QQ plot empirical zTz vs. chi-sq df =14",
  cex = .1)
abline(a = 0, b = 1, col = "blue")
abline(v = quantile(rchisq(100000, df = 14), c(.01, .05, .1, .9, .95, .99)),
       col = c("red", "orange", "skyblue", "skyblue", "orange", "red"))
legend("topleft", legend = c(.9, .95, .99), col = c("skyblue", "orange", "red"),
       lty = 1, title = "percentile")
dev.off()

# QQ plot zTz data vs. MVN01
png("plots/QQ_zTz_data_MVN01.png",
    height = 6, width = 6, units = "in", res = 300)
qqplot(
  x = zTz3_01_mvn,
  y = zTz3,
  main = "QQ plot empirical zTz vs. MVN truncated [0,1]",
  cex = .1)
abline(a = 0, b = 1, col = "blue")
abline(v = quantile(zTz3_01_mvn, c(.01, .05, .1, .9, .95, .99)),
       col = c("red", "orange", "skyblue", "skyblue", "orange", "red"))
legend("topleft", legend = c(.9, .95, .99), col = c("skyblue", "orange", "red"),
       lty = 1, title = "percentile")
dev.off()



# zElev QQ plots


# what is causing the deviation from MVN expectation?
# one thing could be non-normality due to the boundary (or smaller binomial samples)

# this is not the main difference, but even the simulation
# after truncations isn't a great fit:
# I could check for different realized K's after the truncation
zAnc_maize_MVN01 = make_K_calcs(t(mvn_01_maize))
png("plots/K_matrix_diff_MVN01_data_scatterplot.png",
    height = 6, width = 8, units = "in", res = 300)
plot(zAnc_maize$K, zAnc_maize_MVN01$K, cex = 1,
     main = "comparing K matrix before and after truncation [0,1]",
     xlab = "empirical var-cov values going into MVN sim",
     ylab = "var-cov values of sim after truncation",
     xlim = c(0, .05),
     ylim = c(0, .05))
abline(a = 0, b = 1, col = "blue")
dev.off()

# plot a subtraction of the K matrices:
melt(zAnc_maize_MVN01$K - zAnc_maize$K) %>%
  #filter(Var1 != Var2) %>%
  left_join(., meta.pops[, c("ELEVATION", "pop", "LOCALITY")],
            by=c("Var1"="pop")) %>%
  left_join(., meta.pops[, c("ELEVATION", "pop", "LOCALITY")],
            by=c("Var2"="pop")) %>%
  ggplot(data = ., aes(reorder(LOCALITY.x, ELEVATION.x), 
                       y=reorder(LOCALITY.y, ELEVATION.y), 
                       fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1) +  
  ggtitle("K covariance after truncation minus before") +
  xlab('Maize pops low -> high elevation') +
  ylab('Maize pops low -> high elevation')
ggsave("plots/K_matrix_diff_MVN_01_vs_data.png", 
       height = 6, width = 8, 
       units = "in", device = "png")

melt(zAnc_maize$K) %>%
  #filter(Var1 != Var2) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1) +  
  ggtitle("K covariance matrix original")
melt(zAnc_maize_MVN01$K) %>%
  #filter(Var1 != Var2) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(begin = 0, end = 1, direction = 1) +  
  ggtitle("K covariance matrix original after truncation")

# do the alphas change too? some, yes
data.frame(MVN01 = zAnc_maize_MVN01$alpha,
                     empirical = zAnc_maize$alpha,
                     pop = names(zAnc_maize$alpha)) %>%
  left_join(., meta.pops, by = "pop") %>%
  tidyr::gather(., "source", "alpha", c("MVN01", "empirical")) %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = alpha, color = source)) +
  geom_point() +
  ggtitle("Effect of truncating MVN at [0,1] on mean pop mex. ancestry") +
  ylab("Mean mexicana ancestry (alpha)") +
  xlab("maize population low -> high elevation") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/alphas_diff_MVN_01_vs_data.png", 
       height = 6, width = 8, 
       units = "in", device = "png")

# plot zb3_elev against simple elevation
downsample <- c(T, rep(F, 99))
plot(simple_bElev_anc$envWeights[downsample], 
     zb3_elev$zEnv[downsample])
plot(log10(simple_bElev_anc$pval_Env[downsample]), 
     log10(zb3_elev$pval_zEnv[downsample]))
zb3_elev_only <- log10(simple_bElev_anc$pval_Env) < -2 &
  log10(zb3_elev$pval_zEnv) > -1
simple_elev_only <- log10(simple_bElev_anc$pval_Env) > -1 &
  log10(zb3_elev$pval_zEnv) < -2
# plot these:
bind_cols(maize_anc, sites) %>%
  left_join(., zb3_elev[ , c("chr", "pos", "zEnv", "zInt", "pval_zEnv")], by = c("chr", "pos")) %>%
  filter(., zb3_elev_only) %>%
  .[c(T, rep(F, 1000)), ] %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), 
             y = mex_freq, shape = inv4m, 
             group = pos,
             color = zEnv)) +
  geom_point() +
  geom_line() +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset low pvals in zElev but not simple Elev")
ggsave(paste0("plots/a_few_example_loci_zElev_BUT_NOT_simple_elev.png"),
       device = "png",
       width = 10, height = 8, units = "in")

bind_cols(maize_anc, sites, simple_bElev_anc) %>%
  filter(., simple_elev_only) %>%
  .[c(T, rep(F, 50)), ] %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), 
             y = mex_freq, shape = inv4m, 
             group = pos,
             color = envWeights)) +
  geom_point() +
  geom_line() +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset low pvals in simple Elev but not zElev")
ggsave(paste0("plots/a_few_example_loci_simple_elev_BUT_NOT_zElev.png"),
       device = "png",
       width = 10, height = 8, units = "in")

# can I use pvalues then from zEnv??
with(zb3_elev, plot(zEnv[downsample], 
                    pval_zEnv[downsample]))
# much higher correlation in slopes than in pvals
cor(simple_bElev_anc$envWeights, zb3_elev$zEnv)
cor(log10(simple_bElev_anc$pval_Env), 
    log10(zb3_elev$pval_zEnv))
with(zb3_elev, plot(zEnv[downsample], 
                    log10(pval_zEnv[downsample])))
with(simple_bElev_anc, plot(envWeights[downsample], 
                    log10(pval_Env[downsample])))
png("plots/QQ_pval_simple_bElev_vs_uniform.png",
    height = 6, width = 6, units = "in", res = 300)
qqplot(x = runif(10000, 0, 1), 
       y = simple_bElev_anc$pval_Env,
       ylab = "pvalues slope simple bElev",
       main = "QQ plot - pvalues for slope bElev against uniform dist.")
abline(0, 1, col = "blue")
dev.off()
png("plots/QQ_pval_simple_bElev_vs_uniform.png",
    height = 6, width = 6, units = "in", res = 300)
qqplot(x = runif(1000, 0, 1), y = zb3_elev$pval_zEnv,
       ylab = "pvalues slope rotated space zElev",
       main = "QQ plot - pvalues for slope zElev against uniform dist.")
abline(0, 1, col = "blue")
dev.off()
# QQ plot looks reasonable -- but no enrichment for low p-values
table(zb3_elev$pval_zEnv < .01)/nrow(zb3_elev)
table(zb3_elev$pval_zEnv < .001)/nrow(zb3_elev)
table(zb3_elev$pval_zEnv < .0001)/nrow(zb3_elev)
# enrichment for low p's but QQ plot doesn't fit null at all
table(simple_bElev_anc$pval_Env < .01)/nrow(simple_bElev_anc)
table(simple_bElev_anc$pval_Env < .001)/nrow(simple_bElev_anc)
table(simple_bElev_anc$pval_Env < .0001)/nrow(simple_bElev_anc)

# ok let's look at a few high slopes and the models fit to them:
# start with top SNP (in inversion)
max_slope_elev_i <- which.max(simple_bElev_anc$envWeights)
max_slope_zelev_i <- which.max(zb3_elev$zEnv)
# top hit from each model
simple_bElev_anc[max_slope_elev_i,]
zb3_elev[max_slope_elev_i,]
simple_bElev_anc[max_slope_zelev_i,]
zb3_elev[max_slope_zelev_i,]

# ancestry at top outliers - non-rotated:
bind_cols(sites, maize_anc)[c(max_slope_elev_i, max_slope_zelev_i), ] %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = ELEVATION,
             y = mex_freq, 
             shape = inv4m, 
             group = pos)) +
  geom_point(aes(color = LOCALITY)) +
  #geom_line() +
  geom_smooth(method = "lm", color = "black") +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("observed ancestry at top outlier with elevation")

meta.pops.z <- left_join(data.frame(pop = names(zAnc_maize$alpha),
                                    stringsAsFactors = F), 
                          meta.pops, by = "pop")
meta.pops.z$zElevation = zAnc_maize$InvL %*% meta.pops.z$ELEVATION

# also plot in rotated space:
bind_cols(sites, data.frame(zAnc3))[c(max_slope_elev_i, max_slope_zelev_i), ] %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops.z[ , c("pop", "ELEVATION", "LOCALITY", "zElevation")], by = "pop") %>%
  ggplot(aes(x = zElevation,
             y = mex_freq, 
             shape = inv4m, 
             group = pos)) +
  geom_point(aes(color = LOCALITY)) +
  #geom_line() +
  geom_smooth(method = "lm", color = "black") +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("(rotated) ancestry at top outlier with elevation")




# a few example high zTz outliers:
bind_cols(maize_anc, sites) %>%
  mutate(zTz3 = zTz3) %>%
  filter(., zTz3 > quantile(zTz3, .99)) %>%
  .[c(T, rep(F, 250)), ] %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), 
             y = mex_freq, shape = inv4m, 
             group = pos,
             color = zTz3)) +
  geom_point() +
  geom_line() +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset very high zTz3 outliers")
ggsave(paste0("plots/a_few_example_loci_high_zTz_maize.png"),
       device = "png",
       width = 10, height = 8, units = "in")


# try creating a MVBeta (multi-variate beta) using copulas:
# first get draws from a MVN
x = rmvnorm(1000, mean = zAnc_maize$alpha, sigma = zAnc_maize$K)
# then get cumulative distribution for those draws
u = pnorm(x, zAnc_maize$alpha, sqrt(diag(zAnc_maize$K)))
set.seed(1) # set seed then get binomial draws that match the quantiles
# from the simulated normal
k = qbinom(u, 10, zAnc_maize$alpha)
colMeans(k)/10 - zAnc_maize$alpha
#k = qbinom(u, 10, rep(c(.8, .2), 7))
set.seed(1)
k1 = t(sapply(1:1000, function(i) qbinom(u[i, ], 10, zAnc_maize$alpha)))
colMeans(k1)/10 - zAnc_maize$alpha
b = qbeta(u, #thetas estimated in fit_clines.R
          shape1 = zAnc_maize$alpha * coef(m1)[15:28],
          shape2 = (1 - zAnc_maize$alpha) * coef(m1)[1:14])
d = qbeta(u, 
          shape1 = 1, 
          shape2 = 2)
image(cor(x))
image(cor(k)) # matches
image(cor(k1)) # doesn't match -- I'm not sure why. should be the same model?
image(cor(k)-cor(x))
image(cor(b))
image(cor(d))
a <- prob * theta
b <- (1 - prob) * theta
rbeta(n, shape1 = a, shape2 = b)

# try covariance in elevation instead of fitting a linear model
# (because the linear model doesn't constrain the variance to 1, even though we think drift should be)
# zAnc3 saves the zAnc scores for all 14 pops (col), and all loci (rows)
# I want to get one vector zElev (in km elevation):
meta.pops.maize$ELEVATION_km <- meta.pops.maize$ELEVATION/1000
zElev_maize_vector <- c(zAnc_maize$InvL %*% meta.pops.maize$ELEVATION_km)
# center it
zElev_maize_vector_c <- zElev_maize_vector - mean(zElev_maize_vector)
# calculate covariance between observed zAnc and zElev for maize pops
cov_zAnc_zElev_maize <- apply(zAnc3[1:100,], 1, function(row)
  mean((row - mean(row)) * zElev_maize_vector_c))

# try it with centered elevations:
# I want to get one vector zElev:
zElev_maize_vector2 <- zAnc_maize$InvL %*% (meta.pops.maize$ELEVATION_km - mean(meta.pops.maize$ELEVATION_km))
# center it
zElev_maize_vector_c2 <- zElev_maize_vector2 - mean(zElev_maize_vector2)
# calculate covariance between observed zAnc and zElev for maize pops
cov_zAnc_zElev_maize2 <- apply(zAnc3[1:100,], 1, function(row)
  mean((row - mean(row)) * zElev_maize_vector_c2))

# try it with a scalar multiplied by elevation:
# I want to get one vector zElev:
zElev_maize_vector3 <- zAnc_maize$InvL %*% meta.pops.maize$ELEVATION_km/100
# center it
zElev_maize_vector_c3 <- zElev_maize_vector3 - mean(zElev_maize_vector3)
# calculate covariance between observed zAnc and zElev for maize pops
cov_zAnc_zElev_maize3 <- apply(zAnc3[1:100,], 1, function(row)
  mean((row - mean(row)) * zElev_maize_vector_c3))

plot(cov_zAnc_zElev_maize, cov_zAnc_zElev_maize2) # changes relationship to add or subtract from elevation vector
plot(cov_zAnc_zElev_maize, cov_zAnc_zElev_maize3) # invariant to a scalar multiplication

# test expectations/assumptions with small simulation
# the simulation is the pop mean before environmental selection
small_mvn_maize <- mvn_maize[1:1000, ] # just 1000 points
meta.pops.maize$ELEVATION_km_c <- meta.pops.maize$ELEVATION_km - mean(meta.pops.maize$ELEVATION_km)
meta.pops.maize$pop_factor <- 1:14
meta.pops.maize$zElev_c <- zElev_maize_vector_c

# add an allele frequency change that's linear with environment elevation
small_mvn_maize_elev1 <- t(apply(small_mvn_maize, 1, function(row) row + meta.pops.maize$ELEVATION_km_c/5))
# variation with half the slope with environment
small_mvn_maize_elev2 <- t(apply(small_mvn_maize, 1, function(row) row + meta.pops.maize$ELEVATION_km_c/10))
# variation where the zero effect of environment is not constrained to be at the mean observed elevation
small_mvn_maize_elev3 <- t(apply(small_mvn_maize, 1, function(row) row + meta.pops.maize$ELEVATION_km_c/5 + 0.1))
# effect of mean elevation is a slight drop in mexicana ancestry
small_mvn_maize_elev4 <- t(apply(small_mvn_maize, 1, function(row) row + meta.pops.maize$ELEVATION_km_c/5 - 0.1))
# only a mean effect across all samples; no slope with elevation
small_mvn_maize_elev5 <- t(apply(small_mvn_maize, 1, function(row) row + 0.1))


# put neutral and non-neutral sims together in a list
sims_maize_elev <- list(small_mvn_maize, small_mvn_maize_elev1, small_mvn_maize_elev2, 
                        small_mvn_maize_elev3, small_mvn_maize_elev4, small_mvn_maize_elev5)

# calculate zAnc for each simulation:
sims_maize_elev_zAnc <- lapply(sims_maize_elev, function(s)
  t(apply(s, 1, function(l) zAnc(ancFreq = l,
                                 invL = zAnc_maize$InvL,
                                 alpha = zAnc_maize$alpha))))
# plot null model -- should be N(0,1)
hist(sims_maize_elev_zAnc[[1]])
apply(sims_maize_elev_zAnc[[1]], 2, mean) # about 0
apply(sims_maize_elev_zAnc[[1]], 2, var) # about 1

# calculate covariances:
sims_maize_elev_covElev <- lapply(sims_maize_elev_zAnc, function(s)
  apply(s, 1, function(l) mean((l - mean(l)) * zElev_maize_vector_c)))
sd_zElev_maize_vector_c <- sd(zElev_maize_vector_c)
hist(sims_maize_elev_covElev[[1]]/sd_zElev_maize_vector_c)
var(sims_maize_elev_covElev[[1]]/sd_zElev_maize_vector_c)
mean(sims_maize_elev_covElev[[1]]/sd_zElev_maize_vector_c)

lapply(sims_maize_elev_covElev, mean) # 1-3 look good, 4-5 should equal 2 ideally
hist(sims_maize_elev_covElev[[2]]/sd_zElev_maize_vector_c)
hist(sims_maize_elev_covElev[[3]]/sd_zElev_maize_vector_c)
hist(sims_maize_elev_covElev[[4]]/sd_zElev_maize_vector_c)
hist(sims_maize_elev_covElev[[5]]/sd_zElev_maize_vector_c)

# ok, so we don't want to assume the mean effect of elevation in our sample is zero,
# so maybe I can put this back in a linear model and constrain variance on the noise to 1:
# first fit one linear model with all of the generating data
# then try to fit it with limited data (just 14 observations at a locus)
sims_maize_elev_zAnc_d <- lapply(sims_maize_elev_zAnc, function(x)
  x %>%
    data.frame(.) %>%
    data.table::setnames(., meta.pops.maize$pop) %>%
    tidyr::gather(., "pop", "zAnc") %>%
    left_join(., meta.pops.maize[ , c("pop", "pop_factor", "zElev_c")], by = "pop"))
small_mvn_maize_zAnc_d <- sims_maize_elev_zAnc[[1]] %>%
  data.frame(.) %>%
  data.table::setnames(., meta.pops.maize$pop) %>%
  tidyr::gather(., "pop", "zAnc") %>%
  left_join(., meta.pops.maize[ , c("pop", "pop_factor", "zElev_c")], by = "pop")
small_mvn_maize_zAnc_elev1_d <- sims_maize_elev_zAnc[[2]] %>%
  data.frame(.) %>%
  data.table::setnames(., meta.pops.maize$pop) %>%
  tidyr::gather(., "pop", "zAnc") %>%
  left_join(., meta.pops.maize[ , c("pop", "pop_factor", "zElev_c")], by = "pop")


# intercept
zInt = c(zAnc_maize$InvL %*% rep(1, 14))

# run all models:
# this works as expected -- great!
m_elev <- lapply(sims_maize_elev_zAnc_d, function(x)
  map( # quadratic approximation of the posterior MAP
    alist(
      zAnc ~ dnorm(mean = mu, sd = 1),
      mu <- zInt[pop]*b_zInt + zElev_c[pop]*b_zElev,
      b_zElev ~ dnorm(0, 100),
      b_zInt ~ dnorm(0, 10)
    ),
    data = list(zInt = zInt,
                zElev_c = zElev_maize_vector_c,
                zAnc = x$zAnc,
                pop = x$pop_factor)
  )
)
lapply(m_elev, precis)
lapply(m_elev, pairs)

cor(zInt, zElev_maize_vector_c) # 'Neutral' covariance and elevation are highly correlated
# which may mean that a lot of the ancestry covariance is not actually neutral

# can I get estimates for individual loci? I need to label loci in my data
# just look at the distribution of slopes for the neutral simulation:
# make each locus it's own data frame:
small_mvn_maize_zAnc_d <- lapply(1:nrow(sims_maize_elev_zAnc[[1]]), function(i) sims_maize_elev_zAnc[[1]][i, ] %>%
  t(.) %>%
  as.data.frame(.) %>%
  data.table::setnames(., meta.pops.maize$pop) %>%
  tidyr::gather(., "pop", "zAnc") %>%
  left_join(., meta.pops.maize[ , c("pop", "pop_factor", "zElev_c")], by = "pop") %>%
  mutate(snp = i))

m_elev_ind_neutral_loci <- lapply(small_mvn_maize_zAnc_d, function(x)
  map( # quadratic approximation of the posterior MAP
    alist(
      zAnc ~ dnorm(mean = mu, sd = 1),
      mu <- zInt[pop]*b_zInt + zElev_c[pop]*b_zElev,
      b_zElev ~ dnorm(0, 100),
      b_zInt ~ dnorm(0, 10)
    ),
    data = list(zInt = zInt,
                zElev_c = zElev_maize_vector_c,
                zAnc = x$zAnc,
                pop = x$pop_factor)
  )
)
# individual loci are noisy
m_elev_ind_neutral_loci_coef <- do.call(bind_rows, lapply(m_elev_ind_neutral_loci, function(x) x@coef))
precis(m_elev_ind_neutral_loci[[1]])
pairs(m_elev_ind_neutral_loci[[1]])

# a non-neutral scenario:
small_mvn_maize_zAnc_elev1_d <- lapply(1:nrow(sims_maize_elev_zAnc[[2]]), function(i) sims_maize_elev_zAnc[[2]][i, ] %>%
                                   t(.) %>%
                                   as.data.frame(.) %>%
                                   data.table::setnames(., meta.pops.maize$pop) %>%
                                   tidyr::gather(., "pop", "zAnc") %>%
                                   left_join(., meta.pops.maize[ , c("pop", "pop_factor", "zElev_c")], by = "pop") %>%
                                   mutate(snp = i))

m_elev_ind_elev1_loci <- lapply(small_mvn_maize_zAnc_elev1_d, function(x)
  map( # quadratic approximation of the posterior MAP
    alist(
      zAnc ~ dnorm(mean = mu, sd = 1),
      mu <- zInt[pop]*b_zInt + zElev_c[pop]*b_zElev,
      b_zElev ~ dnorm(0, 100),
      b_zInt ~ dnorm(0, 10)
    ),
    data = list(zInt = zInt,
                zElev_c = zElev_maize_vector_c,
                zAnc = x$zAnc,
                pop = x$pop_factor)
  )
)
# individual loci are really too noisy to detect subtle effects like this.
m_elev_ind_elev1_loci_coef <- do.call(bind_rows, lapply(m_elev_ind_elev1_loci, function(x) x@coef))
precis(m_elev_ind_elev1_loci[[1]])
pairs(m_elev_ind_elev1_loci[[1]])
# but are the effects in my data actually on this scale or larger effects?
zAnc3_small <- sample_n(as.data.frame(zAnc3), size = 1000, replace = F) # sample just 1000 loci to get a sense of dist. for outliers
zAnc3_small_d <- lapply(1:nrow(zAnc3_small), function(i) zAnc3_small[i, ] %>%
                                         as.data.frame(.) %>%
                                         data.table::setnames(., meta.pops.maize$pop) %>%
                                         tidyr::gather(., "pop", "zAnc") %>%
                                         left_join(., meta.pops.maize[ , c("pop", "pop_factor", "zElev_c")], by = "pop") %>%
                                         mutate(snp = i))

m_elev_ind_zAnc3_small_loci <- lapply(zAnc3_small_d, function(x)
  map( # quadratic approximation of the posterior MAP
    alist(
      zAnc ~ dnorm(mean = mu, sd = 1),
      mu <- zInt[pop]*b_zInt + zElev_c[pop]*b_zElev,
      b_zElev ~ dnorm(0, 100),
      b_zInt ~ dnorm(0, 10)
    ),
    data = list(zInt = zInt,
                zElev_c = zElev_maize_vector_c,
                zAnc = x$zAnc,
                pop = x$pop_factor)
  )
)
# individual loci are really too noisy to detect subtle effects like this.
m_elev_ind_zAnc3_small_loci_coef <- do.call(bind_rows, lapply(m_elev_ind_zAnc3_small_loci, function(x) x@coef))
precis(m_elev_ind_zAnc3_small_loci[[1]])
pairs(m_elev_ind_zAnc3_small_loci[[1]])
hist(m_elev_ind_neutral_loci_coef$b_zElev, freq = F)
hist(m_elev_ind_zAnc3_small_loci_coef$b_zElev, add = T, freq = F, col = alpha("blue", .3))
hist(m_elev_ind_elev1_loci_coef$b_zElev, add = T, freq = F, col = alpha("green", .3))
summary(m_elev_ind_neutral_loci_coef$b_zElev)
summary(m_elev_ind_zAnc3_small_loci_coef$b_zElev)
summary(m_elev_ind_elev1_loci_coef$b_zElev)
# this set of 1000 loci pretty much look neutral..but then, maybe they are
# what happens with our most extreme outlier?
top_outlier_elev_i <- which(simple_bElev_anc$envWeights == max(simple_bElev_anc$envWeights))
top_outlier_zAnc3_d <- zAnc3[top_outlier_elev_i, ] %>%
                          t(.) %>%
                          as.data.frame(.) %>%
                          tidyr::gather(., "pop", "zAnc") %>%
                          left_join(., meta.pops.maize[ , c("pop", "pop_factor", "zElev_c")], by = "pop") %>%
                          mutate(snp = top_outlier_elev_i)

m_elev_top_outlier_zAnc3 <- map( # quadratic approximation of the posterior MAP
    alist(
      zAnc ~ dnorm(mean = mu, sd = 1),
      mu <- zInt[pop]*b_zInt + zElev_c[pop]*b_zElev,
      b_zElev ~ dnorm(0, 100),
      b_zInt ~ dnorm(0, 10)
    ),
    data = list(zInt = zInt,
                zElev_c = zElev_maize_vector_c,
                zAnc = top_outlier_zAnc3_d$zAnc,
                pop = top_outlier_zAnc3_d$pop_factor)
)
precis(m_elev_top_outlier_zAnc3)
# arg. this should be a very strong positive slope (!), like .6 or .7.
# first check that the populatoin order is the same/as expected in the model.
# then assess other possible problems.

plot(x = meta.pops.maize$ELEVATION[order(meta.pops.maize$ELEVATION)], 
       y = zAnc3[top_outlier_elev_i, ][order(meta.pops.maize$ELEVATION)],
       col = "blue",
       pch = 20)
points(x = meta.pops.maize$ELEVATION[order(meta.pops.maize$ELEVATION)], 
     y = meta.pops.maize$alpha_mex[order(meta.pops.maize$ELEVATION)])

points(x = meta.pops.maize$ELEVATION[order(meta.pops.maize$ELEVATION)], 
       y = maize_anc[top_outlier_elev_i, ][order(meta.pops.maize$ELEVATION)],
       col = "blue")
str(meta.pops.maize)

plot(x = zElev_maize_vector[order(zElev_maize_vector)], 
     y = zAnc3[top_outlier_elev_i, ][order(zElev_maize_vector)],
     col = "blue",
     pch = 20)

# caveat: I haven't truncated simulations to [0,1] scale yet -- another factor influencing sensitivity of the test
plot(x = meta.pops.maize$ELEVATION[order(meta.pops.maize$ELEVATION)], 
       y = maize_anc[top_outlier_elev_i, ][order(meta.pops.maize$ELEVATION)],
       col = "blue")
points(x = meta.pops.maize$ELEVATION[order(meta.pops.maize$ELEVATION)], 
       y = meta.pops.maize$alpha_mex[order(meta.pops.maize$ELEVATION)])
# hmm...this very top outlier just doens't look that impressive in the rotated space..
# partly because it's bounded by [0,1] unrotated, but maybe for other reasons too. I'm sort of at a loss.



# If I mix these 4 simulations across loci, can I estimate a distribution for b_zElev and b_zInt
# I could try to force the mean of b_zInt to be zero if that seems useful

# basic mix of 5 sims allowing them to come from a distribution of b_zElev and b_zInt
mix_sims_maize_elev_zAnc_d <- do.call(rbind, 
                                      lapply(1:length(sims_maize_elev_zAnc_d), function(i)
                                        mutate(sims_maize_elev_zAnc_d[[i]], sim = i)))

# model the mixture of 5 sims:
m_elev_mix <- map( # quadratic approximation of the posterior MAP
  alist(
    zAnc ~ dnorm(mean = mu, sd = 1),
    mu <- zIntercept[pop] * b_zInt[sim] + zElev_c[pop] * b_zElev[sim],
    
    # effect of elevation for mean population
    b_zInt[sim] ~ dnorm(0, 1),
    
    # slope of effect of elevation
    b_zElev[sim] ~ dnorm(0, 10)
  ),
  data = list(zIntercept = zInt,
              zElev_c = zElev_maize_vector_c,
              zAnc = mix_sims_maize_elev_zAnc_d$zAnc,
              pop = mix_sims_maize_elev_zAnc_d$pop_factor,
              sim = mix_sims_maize_elev_zAnc_d$sim)
)
precis(m_elev_mix, depth = 2)
pairs(m_elev_mix)



mix_sims_maize_elev_zAnc_d2 <- do.call(rbind, 
                                       lapply(1:length(sims_maize_elev_zAnc_d), function(i)
                                         mutate(sims_maize_elev_zAnc_d[[i]], sim = i))) %>%
  dplyr::select(-pop) %>%
  left_join(., data.frame(zInt = zInt, 
                          pop_factor = 1:14), by = "pop_factor") %>%
  rename(pop = pop_factor)
m_elev_mix_stan <- map2stan( # quadratic approximation of the posterior MAP
  alist(
    zAnc ~ dnorm(mean = mu, sd = 1),
    mu <- zInt * b_zInt[sim] + zElev_c * b_zElev[sim],
    
    # effect of elevation for mean population
    b_zInt[sim] ~ dnorm(0, 1),
    
    # slope of effect of elevation
    b_zElev[sim] ~ dnorm(0, 10)
  ),
  data = mix_sims_maize_elev_zAnc_d2,
  iter = 1000,
  chains = 2,
  warmup = 200
)
precis(m_elev_mix_stan, depth = 2)
pairs(m_elev_mix_stan)
# plot chains to check for convergence
plot(m_elev_mix_stan, depth = 2)

# model the mixture of 5 sims with stan so it's hierarchical:
m_elev_mix_hier <- map2stan( # quadratic approximation of the posterior MAP
    alist(
      zAnc ~ dnorm(mean = mu, sd = 1),
      mu <- zInt * b_zInt[sim] + zElev_c * b_zElev[sim],
      
      # effect of elevation for mean population
      b_zInt[sim] ~ dnorm(b_zInt_mu, b_zInt_theta),
      b_zInt_mu ~ dnorm(0, 1),
      b_zInt_theta ~ dexp(2),
      
      # slope of effect of elevation
      b_zElev[sim] ~ dnorm(b_zElev_mu, b_zElev_theta),      
      b_zElev_mu ~ dnorm(0, 1),
      b_zElev_theta ~ dexp(2)
    ),
    data = mix_sims_maize_elev_zAnc_d2,
    iter = 1000,
    chains = 2,
    warmup = 500
  )
precis(m_elev_mix_hier)
precis(m_elev_mix_hier, depth = 2)
pairs(m_elev_mix_hier)
# plot chains to check for convergence
plot(m_elev_mix_hier, depth = 2)

# so this works, given many data points.
# How well does it do for just the neutral model if we only give each locus 7 data points?


# and what if I add truncation? Should I just add a logistic after the linear model, and then mvn error (a bit hacky)
# truncate my simulation data. Also maybe make an extreme elevation association example data set.

# compare outliers across populations:
# start with top 1%, or above some mean + 2sd
meta.pops.maize$top_1perc <- apply(maize_anc, 2, function(x) quantile(x, .99))
meta.pops.maize$mean_mex <- apply(maize_anc, 2, mean)
meta.pops.maize$sd2 <- meta.pops.maize$mean_mex + 2*apply(maize_anc, 2, sd)
d_maize_anc_pops <- cbind(sites, maize_anc) %>%
  tidyr::gather(., "pop", "anc", maize_pops) %>%
  left_join(., meta.pops.maize, by = "pop") %>%
  mutate(top1 = anc > top_1perc) %>%
  mutate(top_sd2 = anc > sd2) %>%
  arrange(ELEVATION) %>%
  left_join(., sites_cum, by = c("chr", "pos"))

ind_pop_outliers_by_chr <- lapply(1:10, function(i) 
  #d_maize_anc_pops[c(T, rep(F, 100)),] %>%
  d_maize_anc_pops %>%
  filter(chr == i) %>%
  #mutate(color_by = ifelse(significance == "n.s." & (chr %% 2 == 0), 
  #                         "n.s. - even chrom", significance)) %>%
  ggplot(., aes(pos/10^6, anc, color = top_sd2)) + # plot by position on chromosome (Mb), not relative position
  #ggplot(., aes(pos_cum, anc,
  #              color = color_by, show.legend = F)) +
  geom_hline(data = meta.pops.maize, aes(yintercept = mean_mex), linetype = "dashed", alpha = .5) +
  #geom_point_rast(size = .1) +
  geom_point(size = .1) +
  facet_grid(reorder(LOCALITY, desc(ELEVATION)) ~ .) +
  xlab("position (Mbp)") +
  ylab("mexicana ancestry frequency") +
  ggtitle(paste("Chr", i)) +
  ylim(c(0,1)) +
  #scale_colour_manual(values = c("grey", "skyblue", "orange", "red", "darkgrey")) + 
  scale_colour_manual(values = c("grey", "blue")) + 
  #scale_colour_manual(values = cbPalette[c(1,8)]) + 
  #scale_x_continuous(label = axis_spacing$chr, breaks= axis_spacing$center) +
  theme_classic() +
  coord_cartesian(ylim=c(0,1)) + 
  theme(legend.position = "none") +
  theme(strip.text.y = element_text(angle=0)))
for (i in 1:10){
  ggsave(paste0("plots/ind_pop_outliers_maize_chr", i, ".png"),
         plot = ind_pop_outliers_by_chr[[i]],
         height = 8, width = 16, units = "in",
         device = "png")
  ggsave(paste0("../../hilo_manuscript/figures/ind_pop_outliers_maize_chr", i, ".png"),
         plot = ind_pop_outliers_by_chr[[i]],
         height = 8, width = 16, units = "in",
         device = "png")
}


# what does the null look like from the MVN?
meta.pops.maize$top_1perc_MVN <- apply(mvn_01_maize, 2, function(x) quantile(x, .99))
meta.pops.maize$mean_mex_MVN <- apply(mvn_01_maize, 2, mean)
meta.pops.maize$sd2_MVN <- meta.pops.maize$mean_mex_MVN + 2*apply(mvn_01_maize, 2, sd)
d_mvn_01_maize_pops <- data.frame(mvn_01_maize) %>%
  mutate(pos = 1:nrow(.)) %>% # give arbitrary but unique positions to simulation points
  tidyr::gather(., "pop", "anc", maize_pops) %>%
  left_join(., meta.pops.maize, by = "pop") %>%
  mutate(top1 = anc > top_1perc_MVN) %>%
  mutate(top_sd2 = anc > sd2_MVN) %>%
  arrange(ELEVATION)
# plot SFS of shared outliers across maize populations
d_maize_anc_pops %>%
  group_by(chr, pos) %>%
  summarise(outlier_pops = sum(top_sd2)) %>%
  filter(outlier_pops > 0) %>%
  ggplot(., aes(x = outlier_pops, stat(density))) +
  geom_histogram(fill = "blue", alpha = 0.5) +
  d_mvn_01_maize_pops %>%
  group_by(pos) %>%
  summarise(outlier_pops = sum(top_sd2)) %>%
  filter(outlier_pops > 0) %>%
  #ggplot(., aes(x = outlier_pops, stat(density))) +
  geom_histogram(data = ., aes(x = outlier_pops, stat(density)),
                 fill = "darkgrey", alpha = 0.5) +
  theme_classic()

outliers_maize_bypop <- bind_rows(mutate(d_maize_anc_pops, source = "data"),
                                  mutate(d_mvn_01_maize_pops, source = "mvn_sim") %>%
                                    mutate(inv4m = F)) %>%
  group_by(chr, pos, source, inv4m) %>%
  summarise(outlier_pops = sum(top_sd2))
outliers_maize_bypop %>%
  group_by(source) %>%
  summarise(nonoutlier_perc_loci = sum(outlier_pops == 0)/n(),
            avg_perc_loci_outliers_per_pop = sum(outlier_pops)/n()/14)
outliers_maize_bypop %>%
  filter(outlier_pops > 0) %>%
  ggplot(., aes(x = outlier_pops, stat(density), fill = source)) +
  geom_histogram(position = "dodge2", alpha = .5, binwidth = 0.5) + 
  bind_rows(mutate(d_maize_anc_pops, source = "data") %>%
              filter(!inv4m),
            mutate(d_mvn_01_maize_pops, source = "mvn_sim")) %>%
  group_by(chr, pos, source) %>%
  summarise(outlier_pops = sum(top_sd2)) %>%
  filter(outlier_pops > 0) %>%
  geom_histogram(data = .,
                 position = "dodge2",
                 binwidth = 0.5) +
  scale_fill_manual(values = cbPalette[c(2,1)], name = NULL, labels = c("observed", "simulated")) +
  theme_classic() +
  xlab("Number of populations with high mexicana ancestry") +
  ylab("Density") +
  #theme(legend.position = c(0.8, 0.8)) + # put legend inside plot
  ggtitle("Distribution of shared outliers in maize")
ggsave("plots/histogram_dist_shared_outliers_maize_sd2.png",
       device = "png",
       height = 4, width = 5, units = "in")
ggsave("../../hilo_manuscript/figures/histogram_dist_shared_outliers_maize_sd2.pdf",
       device = "pdf",
       height = 4, width = 5, units = "in")

# plot mexicana too:
d_mexicana_anc_pops <- cbind(sites, mexicana_anc) %>%
  tidyr::gather(., "population", "anc", mexicana_pops) %>%
  left_join(., sites_cum, by = c("chr", "pos"))
d_mexicana_anc_pops %>%
  filter(chr == 4) %>%
  #mutate(color_by = ifelse(significance == "n.s." & (chr %% 2 == 0), 
  #                         "n.s. - even chrom", significance)) %>%
  ggplot(., aes(pos_cum, anc)) + 
  #ggplot(., aes(pos_cum, anc,
  #              color = color_by, show.legend = F)) +
  geom_point(size = .1) +
  facet_grid(population ~ .) +
  xlab("bp position on chromosome") +
  ylab("mexicana ancestry frequency") +
  coord_cartesian(ylim=c(0,1)) + # because xlim limits aren't by default inclusive? awful.
  scale_colour_manual(values = c("grey", "skyblue", "orange", "red", "darkgrey")) + 
  #scale_colour_manual(values = c("grey", "darkgrey")) + 
  #geom_abline(slope = 0, intercept = c(FDRs), linetype = "dashed", alpha = .5, 
  #            color = c("skyblue", "orange", "red")) +
  #scale_x_continuous(label = axis_spacing$chr, breaks= axis_spacing$center) +
  theme(legend.position = "none")
d_mexicana_anc_pops %>%
  filter(population == "pop35") %>%
  filter(chr == 4) %>%
  with(., plot(pos, anc))
d_mexicana_anc_pops %>%
  filter(population == "pop31") %>%
  filter(chr == 4) %>%
  with(., plot(pos_cum, anc))

d_mexicana_anc_pops %>%
  filter(population == "pop31") %>%
  filter(chr == 4) %>%
  summary()
max(d_mexicana_anc_pops$anc)


# TO DO:
# along with the histogram of distribution of shared outliers, look at all combos
# of shared peaks -- do they match with geographical proximity and/or elevation?
# new package to do this: https://github.com/const-ae/ggupset

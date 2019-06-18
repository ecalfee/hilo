# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# TO DO:
# write out generalized least squares comparison to what I'm doing. IF I use AIC it should be AICc for n=14 data points.

rerun_all_models <- F # rerun or just load results of models

# load external scripts
source("../../covAncestry/forqs_sim/k_matrix.R") # import useful functions

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
# separate maize and mexicana data
maize_pops <- unique(meta.pops$pop[meta.pops$zea == "maize"])
mexicana_pops <- unique(meta.pops$pop[meta.pops$zea == "mexicana"])
maize_anc <- all_anc[ , maize_pops]
mexicana_anc <- all_anc[ , mexicana_pops]

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

# calculations
zAnc_mexicana = make_K_calcs(t(mexicana_anc))
zAnc_maize = make_K_calcs(t(maize_anc))

# script to run test cases of ancestry selection models using empirical K matrices
#source("zanc_selection_models.R") # makes plots

# calculate zAnc at each locus. 
zAnc3 <- t(apply(maize_anc, 1, function(l) zAnc(ancFreq = l, 
                                                invL = zAnc_maize$InvL, 
                                                alpha = zAnc_maize$alpha)))

zAnc3_mex <- t(apply(mexicana_anc, 1, function(l) zAnc(ancFreq = l, 
                                                       invL = zAnc_mexicana$InvL, 
                                                       alpha = zAnc_mexicana$alpha)))

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
           envWeights = meta.pops$ELEVATION[meta.pops$zea == "maize"], 
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



# make MVN simulations for maize
n = 100000
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
test_anc <- seq(mean(maize_anc_mean), max(maize_anc_mean), length.out = 10000) # keep in range observed to not divide by 0
test_fdr <- sapply(test_anc, function(x) fdr2(x, data = maize_anc_mean, sims = mvn_01_maize_mean))
png("plots/power_curve_high_mex_in_maize_mvn_sim.png", 
    width = 8, height = 6, units = "in", res = 300)
plot(test_anc, test_fdr, 
     main = "FDR high-mex outliers calculated based on MVN simulation of maize",
     pch = 20, cex = .1)
abline(h=c(.1, .05,.01), col = c("darkgreen", "orange", "red"))
legend("topright", legend = c("FDR = 0.1", "FDR = 0.05", "FDR = 0.01"), col = c("green", "orange", "red"), lty = 1)
dev.off()

# set FDR threshold for low mexicana ancestry in maize
test_anc_low <- seq(min(maize_anc_mean), mean(maize_anc_mean), length.out = 10000) # keep in range observed to not divide by 0
test_fdr_low <- sapply(test_anc_low, function(x) fdr2_low(x, data = maize_anc_mean, sims = mvn_01_maize_mean))
png("plots/power_curve_low_mex_in_maize_mvn_sim.png", 
    width = 8, height = 6, units = "in", res = 300)
plot(test_anc_low, test_fdr_low, # we have no power to detect low mexicana outliers individually 
     main = "FDR low-mex outliers calculated based on MVN simulation of maize",
     pch = 20, cex = .1)
abline(h=c(.1, .05,.01), col = c("darkgreen", "orange", "red"))
legend("topright", legend = c("FDR = 0.1", "FDR = 0.05", "FDR = 0.01"), col = c("green", "orange", "red"), lty = 1)
dev.off()

# calculate FDR thresholds (approx)
FDRs <- sapply(c(0.1, .05, .01), function(p) max(test_anc[test_fdr>p], na.rm = T))


# is seletion for mexicana ancestry a dominant mode of selection?
# what portion of zTz outliers meet 5% FDR threshold for high mex (vs some other mode of selection)



# plot maize ancestry across the genome colored by FDR
d_maize_mean_anc <- sites %>%
  mutate(mean_mexicana_anc_in_maize = maize_anc_mean)
d_maize_mean_anc$significance <- cut(d_maize_mean_anc$mean_mexicana_anc_in_maize,
                                     breaks = c(0, FDRs, 1), include.lowest = F, 
                                     labels = c("n.s.", "FDR < 0.1", "FDR < 0.05", "FDR < 0.01"))
ggplot(d_maize_mean_anc, aes(pos, mean_mexicana_anc_in_maize, color = significance)) +
  geom_point(size = .1) +
  scale_colour_manual(values = c("black", "darkgreen", "orange", "red")) + 
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
simple_bElev_mvn <- apply(mvn_01_maize, 
                       1, function(x)
                       simple_env_regression(ancFreq = x, envWeights = meta.pops.maize$ELEVATION)) %>%
  as.data.frame(t(.))
simple_bElev_anc <- apply(maize_anc, 
                            1, function(x)
                              simple_env_regression(ancFreq = x, envWeights = meta.pops.maize$ELEVATION)) %>%
  as.data.frame(.)
summary(as.data.frame(simple_bElev_mvn)$envWeights)
summary(as.data.frame(simple_bElev_anc)$envWeights)
# get FDR thresholds
test_simple_bElev <- seq(0, max(simple_bElev_anc$envWeights), length.out = 1000) # keep in range observed to not divide by 0
test_fdr_simple_bElev <- sapply(test_simple_bElev, function(x) 
  fdr2(x, data = simple_bElev_anc$envWeights, sims = simple_bElev_mvn$envWeights))
png("plots/power_curve_simple_bElev_in_maize_mvn_sim.png", 
    width = 8, height = 6, units = "in", res = 300)
plot(test_simple_bElev, test_fdr_simple_bElev, 
     main = "FDR pos. slope simple Anc ~ elev calculated based on MVN simulation of maize",
     pch = 20, cex = .1)
abline(h=c(.1, .05,.01), col = c("darkgreen", "orange", "red")) # I think what's happening is only the inversion
# comes out as really enriched .. and it's non-indep so I don't know if that's believable
# ok but overly conservative test because it's correcting for mean mexicana ancestry 
# that is already highly correlated (.84) with elevation
legend("topleft", legend = c("FDR = 0.1", "FDR = 0.05", "FDR = 0.01"), col = c("darkgreen", "orange", "red"), lty = 1)
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
abline(h=c(.1, .05,.01), col = c("darkgreen", "orange", "red")) # I think what's happening is only the inversion
# comes out as really enriched .. and it's non-indep so I don't know if that's believable
# ok but overly conservative test because it's correcting for mean mexicana ancestry 
# that is already highly correlated (.84) with elevation
legend("topleft", legend = c("FDR = 0.1", "FDR = 0.05", "FDR = 0.01"), col = c("darkgreen", "orange", "red"), lty = 1)
dev.off()


FDRs_simple_bElev <- sapply(c(0.1, .05, .01), function(p) max(test_simple_bElev[test_fdr_simple_bElev>p], na.rm = T))
FDRs_simple_bElev_neg <- sapply(c(0.01, .05, .1), function(p) max(test_simple_bElev_neg[test_fdr_simple_bElev_neg<p], na.rm = T))


# plot slopes with environment colored by FDR
d_simple_bElev <- sites %>%
  mutate(slope_elev = simple_bElev_anc$envWeights)
d_simple_bElev$significance <- cut(d_simple_bElev$slope_elev,
                                     breaks = c(-1, FDRs_simple_bElev_neg, FDRs_simple_bElev, 1), include.lowest = F, 
                                     labels = c("FDR < 0.01", "FDR < 0.05", "FDR < 0.1", "n.s.", "FDR < 0.1", "FDR < 0.05", "FDR < 0.01"))

ggplot(d_simple_bElev, 
       aes(pos, slope_elev, color = significance)) +
  geom_point(size = .1) +
  scale_colour_manual(values = c("red", "orange", "darkgreen", "black")) + 
  facet_wrap(~chr, scales = "free_x") +
  geom_abline(slope = 0, intercept = quantile(d_simple_bElev$slope_elev, .99), linetype = "dashed", color = "grey") +
  #geom_abline(slope = 0, intercept = quantile(d_simple_bElev$slope_elev, .01), linetype = "dashed", color = "grey") +
  #geom_abline(slope = 0, intercept = mean(d_simple_bElev$slope_elev), linetype = "dotted", color = "grey")
ggsave("plots/simple_bElev_mex_anc_in_maize_FDR_10_05_01_mvn_sims_whole_genome.png", device = "png",
       height = 6, width = 10, units = "in")

# manhattan style plot
ggplot(d_simple_bElev, 
       aes(pos, slope_elev, color = significance)) +
  geom_point(size = .1) +
  scale_colour_manual(values = c("red", "orange", "darkgreen", "black")) + 
  geom_abline(slope = 0, intercept = quantile(d_simple_bElev$slope_elev, .99), linetype = "dashed", color = "grey") +
  scale_x_continuous( label = axis_spacing$chr, breaks= axis_spacing$center )

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
  summarize(center=( max(pos_cum) + min(pos_cum) ) / 2 )

fdr_bElev = data.frame(FDR = rep(c("FDR < 0.01", "FDR < 0.05", "FDR < 0.1"), 2),
                       cutoff = c(FDRs_simple_bElev, FDRs_simple_bElev_neg),
                       stringsAsFactors = F)
# make wide format                       
d_simple_bElev %>%
  left_join(., sites_cum, by = c("chr", "pos")) %>%
              mutate(even_chr = (chr %% 2 == 0)) %>%
             ggplot(., aes(pos_cum, slope_elev, 
                           color = even_chr, show.legend = F)) +
  geom_point(size = .1) +
  xlab("bp position on chromosome") +
  ylab("slope mexicana ancestry ~ elevation") +
  scale_colour_manual(values = c("grey", "darkgrey")) + 
  geom_abline(slope = 0, intercept = c(FDRs_simple_bElev_neg, FDRs_simple_bElev), linetype = "dashed", color = c("red", "orange", "skyblue", "skyblue", "orange", "red")) +
  geom_abline(slope = 0, intercept = mean(d_simple_bElev$slope_elev), color = "darkgrey", linetype = "dotted") +
  #geom_abline(data = fdr_bElev, aes(intercept = cutoff, slope = 0, color = FDR)) +
  #scale_colour_manual(values = c("grey", "darkgrey", "yellow", "red", "skyblue")) + 
  scale_x_continuous(label = axis_spacing$chr, breaks= axis_spacing$center) +
  theme(legend.position = "none")
ggsave("plots/simple_bElev_mex_anc_in_maize_FDR_10_05_01_mvn_sims_whole_genome_wide.png", device = "png",
       height = 3, width = 12, units = "in")

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


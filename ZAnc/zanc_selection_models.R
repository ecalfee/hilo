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
  dplyr::select(., popN, ELEVATION)
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









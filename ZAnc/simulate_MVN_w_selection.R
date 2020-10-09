# simple simulations of environmental associations + MVN noise
# uses empirical K matrix for either mexicana or maize
# (replaces & extends previous script zEnv_sims.R)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(MASS)
source("ZAnc/k_matrix.R")
source("ZAnc/other_functions.R") # ML_b() calculates maximum likelihood beta estimates
source("ZAnc/lm_env_function.R") # fit simple linear model (regression)

# load variables from Snakefile
# zea = "maize"
# zea = "mexicana"
K_file = snakemake@input[["K"]]
# K_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".K.RData")
meta_file = snakemake@input[["meta_pop"]]
# meta_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pop.meta.RData")
#file_out = snakemake@output[["sim"]]
# file_out = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".MVN.RData")
n_sim = as.integer(snakemake@params[["n_sim"]])
# n_sim = 1000
#fdr_lm = snakemake@output[["fdr_lm"]]
# fdr_lm = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".lmElev.fdr.RData")
Ne = 10000
YESNO = "yes"
plots_prefix = paste0("ZAnc/plots/Ne", Ne, "_", YESNO, "Boot/", zea, "_")

# load data
load(K_file)
load(meta_file) # loads meta_pops
load(fdr_lm) # false discovery rate from empirical data
# calculate elevation in km, mean-centered across pops
meta_pops$elev_km <- meta_pops$ELEVATION/1000
meta_pops$c_elev_km <- meta_pops$elev_km - mean(meta_pops$elev_km)

# functions
# function to minimize neg. loglik and return beta's and loglik
# for environmental selection model
fit_env <- function(y, X, alpha, detK, invK, params){
  o = optim(par = c(0,0),
            fn = function(b) -ll_mvn(y = y,
                                     mu = X %*% b + alpha,
                                     detK = detK,
                                     invK = invK))
  b = broom::tidy(o) %>%
    pivot_wider(data = ., 
                names_from = parameter,
                values_from = value) %>%
    data.table::setnames(params)
  fit = bind_cols(b, broom::glance(o)) %>%
    mutate(ll = -value) %>%
    dplyr::select(-value)
  return(fit)
}

# create some MVN data (noise) that will be
# used by each set of simulations (null & w/ selection)
set.seed(100)
mvn = MASS::mvrnorm(n = n_sim,
                    mu = K$alpha,
                    Sigma = K$K)


# different selection simulations:
# set of slopes (with elevation) and intercepts
# note that the sel[[1]] is the neutral MVN simulation,
# and sel[[2]] and sel[[3]] are a shared selection simulation (with no elevation effect)
# all other simulations have slopes with elevation
slopes = do.call(c, lapply(c(0, 0.1, 0.3, 0.5, 1, 3, -0.1, -0.3, -0.5, -1, -3), rep, 7))
intercepts = rep(c(0, 0.1, 0.3, 0.5, -0.1, -0.3, -0.5), 11)
generating_model = ifelse(intercepts == 0 & slopes == 0, "null",
                          ifelse(slopes == 0, "allpop", "elevation"))
sel_names = 1:length(slopes)

# create set of simulations by adding the simulated
# mvn noise to each expected ancestry frequency,
# under a linear correlation between ancestry and elevation
sel <- lapply(1:length(slopes), function(i)
  t(apply(mvn, 1, function(x) x + intercepts[i] + meta_pops$c_elev_km*slopes[i])))

# truncate simulations to 0-1 ancestry frequency range
# truncate simulated allele freqs within [0,1] range
truncate01 <- function(x){
  t <- x
  t[t < 0] <- 0
  t[t > 1] <- 1
  return(t)
}
sel_trunc <- lapply(sel, truncate01)

# intercept and environment design matrix:
X = cbind(intercept = 1, meta_pops$c_elev_km) %>%
  as.matrix(.)
# solve for inverse of K matrix
invK = solve(K$K)
# solve for determinant of K matrix
detK = det(K$K)

# # get the maximum likelihood beta estimates
# # for each selection simulation
# 
# # how well did we recover simulated values without truncation?
# sel_betas = lapply(sel, function(s)
#   t(apply(s, 1, function(y) 
#     ML_b(y, K$alpha, invK, X))) %>%
#     as.data.frame(.) %>%
#     data.table::setnames(c("b0", "b1")))
# # and with truncation?
# sel_trunc_betas = lapply(sel_trunc, function(s)
#   t(apply(s, 1, function(y) 
#     ML_b(y, K$alpha, invK, X))) %>%
#     as.data.frame(.) %>%
#     data.table::setnames(c("b0", "b1")))

# # finding log likelihood of models
# # sanity check: confirming that optim and analytical solutions are the same:
# # null model
# ll_mvn(y = sel[[1]][1,], 
#        mu = K$alpha, 
#        detK = detK, 
#        invK = invK)
# # selection model
# ll_mvn(y = sel[[1]][1,], 
#        mu = X %*% sel_betas[[1]][1,] + K$alpha, 
#        detK = detK, 
#        invK = invK)
# sel_optims = lapply(sel, function(s) 
#   do.call(rbind, lapply(1:nrow(s), function(i)
#     fit_env(y = s[i,], X = X, alpha = alpha, detK = detK, 
#             invK = invK, params = c("b0", "b1")))))
# sel_trunc_optims = lapply(sel_trunc, function(s) 
#   do.call(rbind, lapply(1:nrow(s), function(i)
#     fit_env(y = s[i,], X = X, alpha = alpha, detK = detK, 
#             invK = invK, params = c("b0", "b1")))))
# 
# # good, very close agreement between optim and analytical optimization to find beta
# plot(sel_optims[[1]]$b0 ~ sel_betas[[1]][ , "b0"])
# plot(sel_trunc_optims[[1]]$b1 ~ sel_trunc_betas[[1]][ , "b1"])
# sapply(1:length(sel_optims), function(i)
#   sum(abs(sel_optims[[i]]$b0 - sel_betas[[i]][ , "b0"]) < .005)
# )
# sapply(1:length(sel_optims), function(i)
#   sum(abs(sel_optims[[i]]$b1 - sel_betas[[i]][ , "b1"]) < .005)
# )
# sapply(1:length(sel_trunc_optims), function(i) # same with truncation
#   sum(abs(sel_trunc_optims[[i]]$b0 - sel_trunc_betas[[i]][ , "b0"]) < .005)
# )
# sapply(1:length(sel_trunc_optims), function(i)
#   sum(abs(sel_trunc_optims[[i]]$b1 - sel_trunc_betas[[i]][ ,"b1"]) < .005)
# )

# merge all sims into one data frame:
# sel_df <- do.call(rbind, lapply(1:length(slopes), function(i)
#   as.data.frame(sel[[i]]) %>%
#     mutate(l = 1:n_sim,
#            intercept = intercepts[i],
#            slope = slopes[i],
#            generating_model = generating_model[i],
#            truncated = F)))
sel_trunc_df <- do.call(rbind, lapply(1:length(slopes), function(i)
  as.data.frame(sel_trunc[[i]]) %>%
    mutate(l = 1:n_sim,
           intercept = intercepts[i],
           slope = slopes[i],
           generating_model = generating_model[i],
           truncated = T)))
#df <- bind_rows(sel_df, sel_trunc_df)
df <- sel_trunc_df

# then fit all selection scenarios (takes a while!):
fits <- c(list(null = fit_null(anc = df[ , names(K$alpha)], 
                                   alpha = K$alpha, invK = invK, detK = detK), 
                   allpop = fit_sel(anc = df[ , names(K$alpha)], 
                                    alpha = K$alpha, invK = invK, detK = detK, 
                                    X = as.matrix(rep(1, length(K$alpha))), 
                                    b_names = "b"), 
                   elevation = fit_sel(anc = df[ , names(K$alpha)], 
                                       alpha = K$alpha, invK = invK, detK = detK, 
                                       X = as.matrix(cbind(intercept = 1, meta_pops$c_elev_km)), 
                                       b_names = c("b0", "b1")))) %>% #,
              # lapply(1:length(K$alpha), function(i)
              #   fit_sel(anc = df[ , names(K$alpha)], 
              #           alpha = K$alpha, invK = invK, detK = detK,
              #           X = as.matrix(1:length(K$alpha) == i), 
              #           b_names = "b"))) %>%
  lapply(., function(x) bind_cols(x, dplyr::select(df, -names(K$alpha))))

#model_names <- c("null", "allpop", "elevation", paste0("onepop", 1:14))
model_names <- c("null", "allpop", "elevation")
names(fits) <- model_names

# using likelihood ratio, compare slope + intercept model (elevation)
# and just intercept model (allpop)
# use Null simulations, not theoretical chisq(df = 1) to set significance
# at fdr for under null at 5% (equivalent to p<.05 2-tailed):
zElev <- left_join(fits[["allpop"]], fits[["elevation"]], 
                   by = c("l", "n", "intercept", "slope", "generating_model", "truncated"), 
                   suffix = c(".all", ".elev")) %>%
  filter(truncated) %>% # only keep truncated data
  dplyr::mutate(diff_ll2 = -2*(ll.all - ll.elev)) # calculate 2*ln(likelihood ratio)

# what happens if I just fit a simple linear model with elevation
# for these simulations?
lmElev <- apply(sel_trunc_df[ , names(K$alpha)], 
                 1, function(x)
                   simple_env_regression(ancFreq = x - K$alpha, # subtract of mean ancestry trend 
                                         envWeights = meta_pops$elev_km)) %>%
  t(.) %>%
  as.data.frame(.) %>%
  bind_cols(dplyr::select(sel_trunc_df, c(l, intercept, slope, generating_model, truncated)), .) %>%
  dplyr::filter(truncated) %>% # only keep truncated
  dplyr::rename(b1 = "envWeights",
               b0 = "(Intercept)") %>%
  dplyr::mutate(abs_slope = abs(b1)) # magnitude of relationship w/ elevation = absolute value of slope

# use only null distribution to determine cutoffs
sig_cutoffs <- bind_rows(
  zElev %>%
  dplyr::filter(slope == 0 & intercept == 0) %>%
  dplyr::summarise(p_0.05 = quantile(diff_ll2, 0.95),
                   p_0.01 = quantile(diff_ll2, 0.99),
                   stat = "zElev"),
  lmElev %>%
    dplyr::filter(slope == 0 & intercept == 0) %>%
    dplyr::summarise(p_0.05 = quantile(abs_slope, 0.95),
                     p_0.01 = quantile(abs_slope, 0.99),
                     stat = "lmElev"))

# count number of significant results per simulated set of slopes and intercepts:
zElev_summary <- zElev %>%
  dplyr::mutate(b1_sign = c("-", "None", "+")[sign(b1) + 2], # sign returns -1, 0, 1 for -/0/+ values
                sig = diff_ll2 > sig_cutoffs$p_0.05[sig_cutoffs$stat == "zElev"]) %>%
  dplyr::group_by(intercept, slope, sig, b1_sign) %>%
  dplyr::summarise(n = n())
lmElev_summary <- lmElev %>%
  dplyr::mutate(b1_sign = c("-", "None", "+")[sign(b1) + 2], # sign returns -1, 0, 1 for -/0/+ values
                sig = abs_slope > sig_cutoffs$p_0.05[sig_cutoffs$stat == "lmElev"]) %>%
  dplyr::group_by(intercept, slope, sig, b1_sign) %>%
  dplyr::summarise(n = n())

# plot how many of the simulations would meet
# 5% FDR cutoff for significance
power_lm <- lmElev_summary %>%
  ggplot(., aes(x = b1_sign, y = n, fill = sig)) + # a few have the wrong inferred sign (oops)
  geom_col() +
  facet_grid(intercept~slope) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab("Sign of inferred slope mexicana anc ~ elevation") +
  scale_fill_manual(values = c("darkgrey", "blue"), labels = c("n.s.", "p < 0.05"), name = "") +
  ggtitle(paste("Power to detect selection w/ elevation using lm with truncation [0,1]", zea))
power_lm
ggsave(paste0(plots_prefix, "power_lm_simulated_MVN_selection.png"),
       plot = power_lm,
       height = 6, width = 8, units = "in", 
       device = "png", dpi = 300)
power_zElev <- zElev_summary %>%
  ggplot(., aes(x = b1_sign, y = n, fill = sig)) + # a few have the wrong inferred sign (oops)
  geom_col() +
  facet_grid(intercept~slope) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab("Sign of inferred slope mexicana anc ~ elevation") +
  scale_fill_manual(values = c("darkgrey", "blue"), labels = c("n.s.", "p < 0.05"), name = "") +
  ggtitle(paste("Power to detect selection w/ elevation using zElev with truncation [0,1]", zea))
power_zElev
ggsave(paste0(plots_prefix, "power_zElev_simulated_MVN_selection.png"),
       plot = power_zElev,
       height = 6, width = 8, units = "in", 
       device = "png", dpi = 300)
# combine into 1 plot to better compare power:
p_true_pos <- bind_rows(zElev_summary %>%
            mutate(stat = "zElev"), 
          lmElev_summary %>%
            mutate(stat = "lmElev")) %>%
  filter(., sig) %>%
  dplyr::mutate(type = ifelse(b1_sign == c("-", "None", "+")[sign(slope) + 2],
    "true positive", # is inferred association in the same direction as simulated slope?
    "false positive")) %>%
  filter(., type == "true positive") %>%
  ggplot(., aes(x = stat, 
                y = n, 
                fill = stat)) +
  geom_col() +
  xlab("") +
  ylab(paste("n sims out of", n_sim)) +
  facet_grid(intercept~slope) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle(paste("True positives p < 0.05 detected from MVN sims + anc ~ elevation [0,1]", zea))
p_true_pos
ggsave(paste0(plots_prefix, "power_true_positives_p0.05_simulated_MVN_selection.png"),
       plot = p_true_pos,
       height = 6, width = 8, units = "in", 
       device = "png", dpi = 300)
p_false_pos <- bind_rows(zElev_summary %>%
            mutate(stat = "zElev"), 
          lmElev_summary %>%
            mutate(stat = "lmElev")) %>%
  filter(., sig) %>%
  dplyr::mutate(type = ifelse(b1_sign == c("-", "None", "+")[sign(slope) + 2],
                              "true positive", # is inferred association in the same direction as simulated slope?
                              "false positive")) %>%
  filter(., type == "false positive") %>%
  ggplot(., aes(x = stat, 
                y = n, 
                fill = stat)) +
  geom_col() +
  xlab("") +
  ylab(paste("n sims out of", n_sim)) +
  facet_grid(intercept~slope) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle(paste("False positives p < 0.05 detected from MVN sims + anc ~ elevation [0,1]", zea))
p_false_pos
ggsave(paste0(plots_prefix, "power_false_positives_p0.05_simulated_MVN_selection.png"),
       plot = p_false_pos,
       height = 6, width = 8, units = "in", 
       device = "png", dpi = 300)

# now what are the signals for these methods along the genome?

#########################################################

bind_rows(zElev_summary %>%
            mutate(stat = "zElev"), 
          lmElev_summary %>%
            mutate(stat = "lmElev")) %>%
  filter(., sig) %>%
  dplyr::mutate(simulated_slope = factor(slope, 
                                            ordered = T,
                                            levels = c(-3, -1, -0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5, 1, 3))) %>%
  ggplot(., aes(x = simulated_slope, 
                y = n, 
                group = stat, 
                color = stat)) +
  geom_point() +
  geom_line() +
  facet_wrap(~intercept)


# make comparable plot for zEnv
# so I only plot sig vs. not sig. for elevation model
# where not sig. includes slope wasn't sig (allpop model wins) and null wins or chisq isn't sig
# so here I always want to give the slope for the elev model
zElev_summary <- best %>%
  dplyr::select(., l, n, intercept, slope, generating_model, slope, best_model, sig, truncated) %>%
  left_join(., dplyr::select(fits[["elevation"]], l, n, intercept, slope, truncated, b0, b1),
            by = c("l", "n", "intercept", "slope", "truncated")) %>%
  dplyr::mutate(b1_sign = ifelse(is.na(b1), "None",  # what sign does the inferred slope have?
                                 ifelse(b1 > 0, "+", "-")),
                b0_sign = ifelse(is.na(b0), "None", 
                                 ifelse(b0 > 0, "+", "-")),
                sig = (sig & best_model == "elevation")) %>%
  dplyr::group_by(intercept, slope, truncated, sig, b1_sign, b0_sign) %>%
  dplyr::summarise(n = n())

power_zElev <- zElev_summary %>%
  filter(truncated) %>%
  ggplot(., aes(x = b1_sign, y = n, fill = sig)) + # a few have the wrong inferred sign (oops)
  geom_col() +
  facet_grid(intercept~slope) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
  ggtitle(paste("Power to detect selection w/ elevation using ll-ratio with truncation [0,1]", zea))
ggsave(paste0(plots_prefix, "power_zElev_simulated_MVN_b1_b0.png"),
       plot = power_zElev,
       height = 6, width = 8, units = "in", 
       device = "png", dpi = 300)

# conclusion: zElev model has higher power than simple slope model, especially
# for changes of 0.5-1 and -0.3 to -0.5 proportion mexicana over 1km elevation gain
# smaller slopes both methods have very low power
# and larger slopes they both do very well
# having on avg. higher or lower intercept (mexicana %)
# makes a big diff in power to detect a slope (likely b/c of bounds [0,1])
# and 'significant' slopes that are very small under the loglik-ratio test are unreliable (can be false + in wrong direction if mexicana ancestry is elevated and abs(slope) <= 0.1)









# what % are 'sig' under an elevation model?
d <- do.call(bind_rows, lapply(1:3, function(i) fits[[i]] %>%
                                 mutate(fit_model = names(fits)[i]))) %>%
  pivot_wider(data =., id_cols = c(l, n, intercept, slope, generating_model, truncated),
              names_from = fit_model,
              values_from = c("RSS", "ll", "AICc", "k", "b0", "b1"))
# find best model (by AICc). Also test if sig. better than null.
best <- do.call(bind_rows, lapply(1:3, function(i) fits[[i]] %>%
                                 mutate(fit_model = names(fits)[i]))) %>%
  dplyr::group_by(l, n, intercept, slope, generating_model, truncated) %>%
  dplyr::summarise(best_AICc = min(AICc),
                   best_model = fit_model[which.min(AICc)],
                   best_ll = ll[which.min(AICc)],
                   best_k = k[which.min(AICc)],
                   b0 = b0[which.min(AICc)],
                   b1 = b1[which.min(AICc)]) %>%
  left_join(., 
            dplyr::select(fits[["null"]], l, n, slope, intercept, generating_model, truncated, ll, AICc), 
            by = c("l", "n", "slope", "intercept", "generating_model", "truncated")) %>%
  dplyr::mutate(chisq_critial_value01 = qchisq(0.99, df = best_k), # calculate critical value for chisq  significance at alpha = 0.01
                ll_ratio = -2*(ll - best_ll),
                sig = ll_ratio > chisq_critial_value01)
  
# how many of the simulated values are significant for diff. b0 and b1?
best_summary <- best %>%
  dplyr::mutate(b1_sign = ifelse(is.na(b1), "None",  # what sign does the inferred slope have?
                ifelse(b1 > 0, "+", "-")),
                b0_sign = ifelse(is.na(b0), "None", 
                                 ifelse(b0 > 0, "+", "-"))) %>%
  dplyr::group_by(intercept, slope, truncated, best_model, sig, b1_sign, b0_sign) %>%
  dplyr::summarise(n = n())


# using chisq log likelihood ratio test
# how frequently would we detect a sig. outlier
# for different true slopes b1 and intercept b0?
# using the truncated simulated data [0,1]
power_best_model <- best_summary %>%
  filter(truncated) %>%
  ggplot(., aes(x = best_model, y = n, fill = sig)) + # a few have the wrong inferred sign (oops)
  geom_col() +
  facet_grid(intercept~slope) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle(paste("Power to detect selection w/ ll-ratio test with truncation [0,1]", zea))
ggsave(paste0(plots_prefix, "power_best_model_simulated_MVN_b1_b0.png"),
       plot = power_best_model,
       height = 6, width = 8, units = "in", 
       device = "png", dpi = 300)


# how strong is the winner's curse here?
# what are the estimated b1's for the loci where we detect a slope?
#filter(zElev_summary, slope == 0.5) %>% ggplot(., aes(x = b1)) + geom_histogram() + facet_wrap(~sig)
plot_sig_slopes_zElev <- best %>%
  dplyr::select(., l, n, intercept, slope, generating_model, slope, best_model, sig, truncated) %>%
  left_join(., dplyr::select(fits[["elevation"]], l, n, intercept, slope, truncated, b0, b1),
            by = c("l", "n", "intercept", "slope", "truncated")) %>%
  dplyr::mutate(b1_sign = ifelse(is.na(b1), "None",  # what sign does the inferred slope have?
                                 ifelse(b1 > 0, "+", "-")),
                b0_sign = ifelse(is.na(b0), "None", 
                                 ifelse(b0 > 0, "+", "-")),
                sig = (sig & best_model == "elevation")) %>%
  filter(slope >= 0 & truncated) %>%
  ggplot(., aes(x = b1, fill = sig)) +
  geom_histogram() +
  facet_grid(intercept~slope) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle(paste("Model zElev estimated slopes from MVN sim truncated [0,1]", zea))
plot_sig_slopes_zElev
ggsave(paste0(plots_prefix, "slopes_zElev_simulated_MVN_b1_b0.png"),
       plot = plot_sig_slopes_zElev,
       height = 6, width = 8, units = "in", 
       device = "png", dpi = 300)


# make equivalent plot for linear model estimates:
plot_sig_slopes_lm <- fits_lm %>%
  filter(slope >= 0 & truncated) %>%
  ggplot(., aes(x = envWeights, fill = sig)) +
  geom_histogram() +
  facet_grid(intercept~slope) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ggtitle(paste("Model lm() estimated slopes from MVN sim truncated [0,1]", zea))
# geom_vline(data = fits_lm %>%
#                 filter(slope >= 0 & truncated) %>%
#                 group_by(slope, sig, intercept) %>%
#                 summarise(envWeights = mean(envWeights)),
#               aes(xintercept = envWeights, color = sig))
plot_sig_slopes_lm
ggsave(paste0(plots_prefix, "slopes_lm_simulated_MVN_b1_b0.png"),
       plot = plot_sig_slopes_lm,
       height = 6, width = 8, units = "in", 
       device = "png", dpi = 300)

# plot just a couple of the simulated ancestry frequencies
df %>%
  filter(slope == 0, intercept == 0.3, truncated == T) %>%
  sample_n(9) %>%
  pivot_longer(data = ., cols = starts_with("pop"), values_to = "anc", names_to = "pop") %>%
  left_join(meta_pops, by = "pop") %>%
  ggplot(., aes(x = ELEVATION, y = anc, color = l, group = l)) +
  geom_line() +
  facet_wrap(~l) +
  ylim(0:1)

# how many sig. sites in the real data?
# are you convinced by the model?
# e.g. inv4m?
# fit just 1 site in inversion:


# goal: to compare power to detect associations with environment
# using zEnv or simple linear models

# So lower AICc is better. How many true elevation models 
# w/ each choice of parameters have 
# AICc elevation model < AICc null (and < AICc other models? esp. intercept only model)
# Q. is RSS under the null = ZAnc statistic? Check on this.
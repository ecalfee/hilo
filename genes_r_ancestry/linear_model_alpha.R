# in this script I fit a linear model for the admixture proportion alpha
# using estimated allele frequencies for allopatric and sympatric populations
library(dplyr)
library(bedr)
library(rethinking)
# use multiple cores for model fits
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# color blind palette:
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# maize vs. mexicana (sympatric) vs. parviglumis
#scale_color_manual(values = cbPalette[c(2,4,8)])
# maize vs. mexicana (allopatric)
#scale_color_manual(values = cbPalette[c(7,6)])


pops = c("allo.maize", "allo.mexicana", "symp.maize", "symp.mexicana", "parv")
pops_extra = c("pop35", "trip")
# paths to input allele frequencies
path_allele_freqs = "../within_ancestry_diversity/results/allele_freq/pass2_alloMAIZE"

# get sites SNP data
sites_file = paste0("../variant_sites/results/pass2_alloMAIZE/all_regions_every_100th.var.sites")

sites0 <- read.table(sites_file, 
                            stringsAsFactors = F, header = F) %>%
  data.table::setnames(c("chr", "pos", "major", "minor"))
# add recombination rates
# 1cM resolution:
map_pos_1cM <- read.table("results/map_pos_1cM_windows.txt",
            header = T, sep = "\t") %>%
  mutate(chr = as.character(chr))

sites0_bed <- sites0 %>%
  mutate(start = pos - 1, end = pos) %>%
  dplyr::select(c("chr", "start", "end")) %>%
  mutate(chr = as.character(chr))
  
sitesR <- bedr(
  engine = "bedtools", 
  input = list(a = sites0_bed, b = map_pos_1cM), 
  method = "map", 
  params = "-c 5,6,4 -o collapse -g ../data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths",
  check.chr = F
) %>%
  data.table::setnames(c("chr", "start", "end", "mean_cM_Mb", "bin_r5", "window")) %>%
  mutate(mean_cM_Mb = as.numeric(mean_cM_Mb))


allele_freq <- do.call(cbind,
                       lapply(pops, function(POP)
  read.table(paste0(path_allele_freqs, "/", POP, "/all_regions_every_100th.freqs.txt"),
             stringsAsFactors = F, header = F))) %>%
    data.table::setnames(pops) %>%
  cbind(sites0, sitesR[ , c("mean_cM_Mb", "bin_r5", "window")], .)

allele_freq_extra <- do.call(cbind,
                       lapply(pops_extra, function(POP)
                         read.table(paste0(path_allele_freqs, "/", POP, "/all_regions_every_100th.freqs.txt"),
                                    stringsAsFactors = F, header = F))) %>%
  data.table::setnames(pops_extra) %>%
  cbind(allele_freq, .)

f <- allele_freq_extra[complete.cases(allele_freq), ] %>%
  mutate(bin_r5 = as.factor(bin_r5)) %>%
  mutate(bin_r5_value = as.numeric(bin_r5))
# simple test look good
with(f, lm(symp.maize ~ allo.maize + allo.mexicana)) # mostly maize
with(f, lm(symp.mexicana ~ allo.maize + allo.mexicana)) # mostly mex
with(f, lm(parv ~ allo.maize + allo.mexicana)) # about half and half

# pop35 is smaller but may be a less admixed proxy for mexicana ancestry
with(f, lm(allo.mexicana ~ allo.maize + pop35))
with(f, lm(pop35 ~ allo.maize + allo.mexicana))

m0_maize <- map( # quadratic approximation of the posterior MAP
  alist(
    symp.maize ~ dnorm(p, theta),
    #symp.maize ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    p <- alpha*allo.mexicana + (1 - alpha)*allo.maize,
    alpha ~ dunif(0, 1),
    theta ~ dunif(0, 10)
  ),
  data = f[(f$allo.maize != f$allo.mexicana), ],
  start = list(theta = 2, alpha = 0.3))
precis(m0_maize)
pairs(m0_maize)


m1_maize <- map( # quadratic approximation of the posterior MAP
  alist(
    symp.maize ~ dnorm(p, theta),
    #symp.maize ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    p <- logistic(alpha0 + mean_cM_Mb*b_r)*allo.mexicana + (1 - logistic(alpha0 + mean_cM_Mb*b_r))*allo.maize,
    #logit(alpha) <- alpha0 + mean_cM_Mb*b_r,
    alpha0 ~ dnorm(0, 5),
    b_r ~ dnorm(0, 5),
    theta ~ dunif(0, 10)
  ),
  data = f[(f$allo.maize != f$allo.mexicana), ],
  start = list(theta = 2, alpha0 = 0, b_r = 0))
precis(m1_maize)
pairs(m1_maize)
logistic(coef(m1_maize)["alpha0"]+coef(m1_maize)["b_r"]*c(.01, .5, 3, 100))

# pretty slow but I can alternatively formulate this same model in stan
m1_maize_stan <- map2stan( # quadratic approximation of the posterior MAP
  alist(
    symp.maize ~ dnorm(p, theta),
    #symp.maize ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    p <- alpha*allo.mexicana + (1 - alpha)*allo.maize,
    logit(alpha) <- alpha0 + mean_cM_Mb*b_r,
    alpha0 ~ dnorm(0, 5),
    b_r ~ dnorm(0, 5),
    theta ~ dunif(0, 10)
  ),
  data = f[(f$allo.maize != f$allo.mexicana), ],
  start = list(theta = 2, alpha0 = 0, b_r = 0),
  chains = 2, iter = 1000, warmup = 500)
precis(m1_maize_stan)
pairs(m1_maize_stan)
logistic(coef(m1_maize_stan)["alpha0"]+coef(m1_maize_stan)["b_r"]*c(.01, .5, 3, 100))


# by bin, not r directly.
# errors running (!)
m1_maize_rbin <- map( # quadratic approximation of the posterior MAP
  alist(
    symp.maize ~ dnorm(p, theta),
    #symp.maize ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    #p <- alpha*allo.mexicana + (1 - alpha)*allo.maize,
    #logit(alpha) <- alpha0 + bin_r5_value*b_rbin,
    p <- logistic(alpha0 + b_rbin[bin_r5_value])*allo.mexicana + 
      (1 - logistic(alpha0 + b_rbin[bin_r5_value]))*allo.maize,
    alpha0 ~ dnorm(0, 5),
    # each bin has it's own fixed effect (no imposed hierarchy)
    b_rbin[bin_r5_value] ~ dnorm(0, 5),
    theta ~ dunif(0, 10)
  ),
  data = f[(f$allo.maize != f$allo.mexicana), ],
  start = list(theta = 2, alpha0 = 0, b_rbin = 0))
precis(m1_maize_rbin, depth = 2)
pairs(m1_maize_rbin, depth = 2)
logistic(coef(m1_maize_rbin)["alpha0"] + coef(m1_maize_rbin)[3:7])
# same as m0_maize_rbin below (easy check)

# run each bin separately:
m0_maize_rbin <- lapply(1:5, function(i) 
                        map( # quadratic approximation of the posterior MAP
  alist(
    symp.maize ~ dnorm(p, theta),
    #symp.maize ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    p <- alpha*allo.mexicana + (1 - alpha)*allo.maize,
    alpha ~ dunif(0, 1),
    theta ~ dunif(0, 10)
  ),
  data = f[f$bin_r5_value == i & (f$allo.maize != f$allo.mexicana), ],
  start = list(theta = 2, alpha = 0.3)))
lapply(m0_maize_rbin, precis)

# mexicana:
m0_mexicana <- map( # quadratic approximation of the posterior MAP
  alist(
    symp.mexicana ~ dnorm(p, theta),
    p <- alpha*allo.mexicana + (1 - alpha)*allo.maize,
    alpha ~ dunif(0, 1),
    theta ~ dunif(0, 10)
  ),
  data = f[(f$allo.maize != f$allo.mexicana), ],
  start = list(theta = 2, alpha = 0.3))
precis(m0_mexicana)
pairs(m0_mexicana)


m1_mexicana <- map( # quadratic approximation of the posterior MAP
  alist(
    symp.mexicana ~ dnorm(p, theta),
    p <- logistic(alpha0 + mean_cM_Mb*b_r)*allo.mexicana + (1 - logistic(alpha0 + mean_cM_Mb*b_r))*allo.maize,
    alpha0 ~ dnorm(0, 5),
    b_r ~ dnorm(0, 5),
    theta ~ dunif(0, 10)
  ),
  data = f[(f$allo.maize != f$allo.mexicana), ],
  start = list(theta = 2, alpha0 = 0)) #, b_r = 0))
precis(m1_mexicana)
pairs(m1_mexicana)
logistic(coef(m1_mexicana)["alpha0"]+coef(m1_mexicana)["b_r"]*c(.01, .5, 3, 100))

# by bin, not r directly.
# errors running (!)
m1_mexicana_rbin <- map( # quadratic approximation of the posterior MAP
  alist(
    symp.mexicana ~ dnorm(p, theta),
    p <- logistic(alpha0 + b_rbin[bin_r5_value])*allo.mexicana + 
      (1 - logistic(alpha0 + b_rbin[bin_r5_value]))*allo.maize,
    alpha0 ~ dnorm(0, 5),
    # each bin has it's own fixed effect (no imposed hierarchy)
    b_rbin[bin_r5_value] ~ dnorm(0, 5),
    theta ~ dunif(0, 10)
  ),
  data = f[(f$allo.maize != f$allo.mexicana), ],
  start = list(theta = 2, alpha0 = 0, b_rbin = 0))
precis(m1_mexicana_rbin, depth = 2)
pairs(m1_mexicana_rbin, depth = 2)
logistic(coef(m1_mexicana_rbin)["alpha0"] + coef(m1_mexicana_rbin)[3:7])

# I don't expect to see any relationship between ancestry and r in parviglumis
m1_parv <- map( # quadratic approximation of the posterior MAP
  alist(
    parv ~ dnorm(p, theta),
    p <- logistic(alpha0 + mean_cM_Mb*b_r)*allo.mexicana + (1 - logistic(alpha0 + mean_cM_Mb*b_r))*allo.maize,
    alpha0 ~ dnorm(0, 5),
    b_r ~ dnorm(0, 5),
    theta ~ dunif(0, 10)
  ),
  data = f[(f$allo.maize != f$allo.mexicana), ],
  start = list(theta = 2, alpha0 = 0)) #, b_r = 0))
precis(m1_parv)
pairs(m1_parv)
logistic(coef(m1_parv)["alpha0"]+coef(m1_parv)["b_r"]*c(.01, .5, 3, 100))

m1_parv_rbin <- map( # quadratic approximation of the posterior MAP
  alist(
    parv ~ dnorm(p, theta),
    p <- logistic(alpha0 + b_rbin[bin_r5_value])*allo.mexicana + 
      (1 - logistic(alpha0 + b_rbin[bin_r5_value]))*allo.maize,
    alpha0 ~ dnorm(0, 5),
    # each bin has it's own fixed effect (no imposed hierarchy)
    b_rbin[bin_r5_value] ~ dnorm(0, 5),
    theta ~ dunif(0, 10)
  ),
  data = f[(f$allo.maize != f$allo.mexicana), ],
  start = list(theta = 2, alpha0 = 0, b_rbin = 0))
precis(m1_parv_rbin, depth = 2)
pairs(m1_parv_rbin, depth = 2)
logistic(coef(m1_parv_rbin)["alpha0"] + coef(m1_parv_rbin)[3:7])

# use pop35 instead of allopatric mexicana group:
m0_maize_35 <- map( # quadratic approximation of the posterior MAP
  alist(
    symp.maize ~ dnorm(p, theta),
    #symp.maize ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    p <- alpha*pop35 + (1 - alpha)*allo.maize,
    alpha ~ dunif(0, 1),
    theta ~ dunif(0, 10)
  ),
  data = f[(!is.na(f$pop35) & f$allo.maize != f$pop35), ],
  start = list(theta = 2, alpha = 0.5))
precis(m0_maize_35)
pairs(m0_maize_35)

m1_maize_35 <- map( # quadratic approximation of the posterior MAP
  alist(
    symp.maize ~ dnorm(p, theta),
    #symp.maize ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    p <- logistic(alpha0 + mean_cM_Mb*b_r)*pop35 + (1 - logistic(alpha0 + mean_cM_Mb*b_r))*allo.maize,
    #logit(alpha) <- alpha0 + mean_cM_Mb*b_r,
    alpha0 ~ dnorm(0, 5),
    b_r ~ dnorm(0, 5),
    theta ~ dunif(0, 10)
  ),
  data = f[(!is.na(f$pop35) & f$allo.maize != f$pop35), ],
  start = list(theta = 2, alpha0 = 0, b_r = 0))
precis(m1_maize_35)
pairs(m1_maize_35)
logistic(coef(m1_maize_35)["alpha0"]+coef(m1_maize_35)["b_r"]*c(.01, .5, 3, 100))

# model comparison
rethinking::compare(m0_maize, m1_maize, m1_maize_rbin)


# make plots of model results:
bin_labels <- f %>%
  group_by(bin_r5) %>%
  summarise(r = mean(mean_cM_Mb))
png("plots/linear_model_alpha_estimate_by_r.png",
    height = 6, width = 8, units = "in", res = 300)
curve(logistic(coef(m1_maize)["alpha0"] + coef(m1_maize)["b_r"]*x), 
      from = 0, to = 100,
      n = 100,
      xlim = c(0,100),
      ylim = c(0,1),
      ylab = "mexicana ancestry proportion",
      xlab = "recombination rate cM/Mb at 1cM resolution",
      main = "linear model fit: effect of recombination on alpha",
      col = cbPalette[2])
curve(logistic(coef(m1_mexicana)["alpha0"] + coef(m1_mexicana)["b_r"]*x), 
      from = 0, to = 100,
      col = cbPalette[4],
      add = T)
curve(logistic(coef(m1_parv)["alpha0"] + coef(m1_parv)["b_r"]*x), 
      from = 0, to = 100,
      col = cbPalette[8],
      add = T)
points(bin_labels$r, logistic(coef(m1_mexicana_rbin)["alpha0"] + coef(m1_mexicana_rbin)[3:7]),
       add = T, col = cbPalette[4], pch = c(0,1,2,5,6))
points(bin_labels$r, logistic(coef(m1_maize_rbin)["alpha0"] + coef(m1_maize_rbin)[3:7]),
       add = T, col = cbPalette[2], pch = c(0,1,2,5,6))
points(bin_labels$r, logistic(coef(m1_parv_rbin)["alpha0"] + coef(m1_parv_rbin)[3:7]),
       add = T, col = cbPalette[8], pch = c(0,1,2,5,6))
legend("right", c("sympatric maize", "sympatric mexicana", "parviglumis", levels(bin_labels$bin_r5)), 
       pch = c(16, 16, 16, 0, 1, 2, 5, 6), 
       col = c(cbPalette[c(2,4,8)], rep("black", 5)),
       cex = .5)
dev.off()       

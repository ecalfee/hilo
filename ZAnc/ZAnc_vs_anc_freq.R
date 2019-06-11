library(scales)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(MASS)
library(crch) # for censored normal dist.

source("../../covAncestry/forqs_sim/k_matrix.R") # import useful functions

d <- read.table("results/pass1_maize/pop.ZAnc", stringsAsFactors = F,
                header = T)
summary(d)
quantile(d$mean_anc, c(.5, .95, .99, .999))
quantile(d$ZAnc, c(.5, .95, .99, .999))
png(file = "plots/corr_ZAnc_and_mean_ancestry.png",
    width = 10, height = 10, units = "in", res = 300)
with(d, plot(mean_anc, ZAnc, 
             cex = 1,
             col = alpha(ifelse(mean_anc > 0.55 & ZAnc > 13, "red", "blue"), .05),
             pch = 18,
             main = "99% ZAnc and Mean Ancestry across pops"))
abline(h = mean(d$ZAnc), v = mean(d$mean_anc))
dev.off()


mean(d$mean_anc)
sd(d$mean_anc)

# to do: what does the ancestry cov. matrix look like?
# for weird low ZAnc but moderate mean_anc (0.17 ish)

# plot what individual populations are doing
# first get individual population meta data
pass1_allo4Low <- read.table("../data/pass1_allo4Low_ids.txt", stringsAsFactors = F, 
                             header = T, sep = "\t")

# combine sample meta data
# add elevation
meta <- read.table("../data/riplasm/gps_and_elevation_for_sample_sites.txt",
                   stringsAsFactors = F, header = T, sep = "\t") %>%
  dplyr::select(., popN, ELEVATION) %>%
  left_join(pass1_allo4Low, ., by = "popN") %>%
  mutate(., group = paste(symp_allo, zea, sep = "_"))



dir_anc <- "../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/output_noBoot/anc"
filter(meta, group == "sympatric_maize" & est_coverage >= .05) %>%
  group_by(popN)%>%
  #summarise(., count = n())
  arrange(popN) %>%
  unique(.) %>%
  filter(!(popN %in% c(361, 369, 371)))
  #do.call(cbind,
pops.symp.maize <- c(360, 362, 363, 365:368, 370, 372:374)
pops.low.symp.maize <- c(361, 369, 371)
a <- do.call(cbind,
  lapply(pops.symp.maize, function(i)
         read.table(paste0(dir_anc, "/pop", i, ".anc.freq"))))
colnames(a) <- paste0("pop", pops.symp.maize)
summary(a)
# add in very low mexicana ancestry pops:
a_low <- do.call(cbind,
                 lapply(pops.low.symp.maize, function(i)
                 read.table(paste0("../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/output_high_maize_prior_noBoot/anc",
                                   "/pop", i, ".anc.freq"))))
colnames(a_low) <- paste0("pop", pops.low.symp.maize)
summary(a_low)
a2 <- cbind(a, a_low)
# now make more meaningful population labels
counts <- meta %>%
  filter(est_coverage >= 0.05) %>%
  group_by(popN, group, LOCALITY, ELEVATION, zea, symp_allo) %>%
  summarize(count = n())
a2.pops <- data.frame(popN = as.integer(substr(colnames(a2), 4, 6)),
                      alpha = colMeans(a2),
                      var_observed = apply(a2, 2, var)) %>%
  left_join(., counts, by = "popN") %>%
  mutate(var_npq = count*2*alpha*(1-alpha)/(count*2)^2) %>%
  mutate(var_obs_corrected = var_observed*(count*2-1)/(count*2)) # I don't think needed here
colnames(a2) <- a2.pops$LOCALITY

ggplot(a2.pops, aes(x = ELEVATION, y = alpha, color = LOCALITY)) +
  geom_point() +
  ggtitle("Mean of local ancestry calls per pop by elevation")
ggsave("plots/mean_anc_from_hmm_by_elevation.png",
       height = 8, width = 10, device = "png", units = "in")

# what is the theoretical minimum of ZAnc based on the ancestry variance-covariance matrix?
# calculate variance-covariance matrix
a2.alpha <- apply(a2, 2, mean) # mean per pop
a2.cov <- cov(a2)
# calc min ZAnc assuming all pops have zero ancestry
ZAnc(ancFreq = rep(0, length(a2.alpha)), invL = calcInvL(a2.cov), alpha = a2.alpha)
a2.K <- calcK(t(a2), alpha = a2.alpha)
ZAnc(ancFreq = rep(0, length(a2.alpha)), invL = calcInvL(a2.K), alpha = a2.alpha)
# ok why do I get ZAnc's that are lower than the theoretical min?
a2.ZAnc <- apply(a2, 1, function(x) 
  ZAnc(ancFreq = x, invL = calcInvL(a2.K), alpha = a2.alpha))
summary(a2.ZAnc) # how am I getting ZAnc below what I should theoretically get for zero ancestry?
length(which(a2.ZAnc < -3.6))/length(a2.ZAnc) # nearly 10% ZAnc stats are too low
a2[(which(a2.ZAnc == min(a2.ZAnc))),]
a2[(which(a2.ZAnc == min(a2.ZAnc))),] - a2.alpha

calcInvL(a2.K) %*% t(a2[(which(a2.ZAnc == min(a2.ZAnc))), ] - a2.alpha)
sum(calcInvL(a2.K) %*% t(a2[(which(a2.ZAnc == min(a2.ZAnc))), ] - a2.alpha))
ZAnc(ancFreq = t(a2[(which(a2.ZAnc == min(a2.ZAnc))), ]), invL = calcInvL(a2.K), alpha = a2.alpha)
zAnc(ancFreq = t(a2[(which(a2.ZAnc == min(a2.ZAnc))), ]), invL = calcInvL(a2.K), alpha = a2.alpha)

# look at zAnc, i.e. before summing
a2.zAnc <- do.call(cbind,
                   lapply(1:nrow(a2), function(i) 
  zAnc(ancFreq = t(a2[i, ]), invL = calcInvL(a2.K), alpha = a2.alpha)))
# did the rotation work as we'd expect? Looks like yes
a2.zAnc.cov <- cov(t(a2.zAnc))
# plot variances and covariances:
melt(a2.zAnc.cov) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("var-cov of maize pop ancestries after zAnc transformation")
ggsave("plots/variance_covariance_matrix_zAnc_transformed_maize.png", 
       width = 8, height = 8, units = "in",
       device = "png")
png("plots/ZAnc_maize_hist.png")
hist(apply(a2.zAnc, 2, sum), main = "histogram of ZAnc score - maize", freq = F)
curve(expr = dnorm(x, 0, sqrt(14)), from = -25, to = 25, add = T, col = "darkblue")
dev.off()

# what do these very low ZAnc score loci look like?
a2.lowZ <- a2[which(a2.ZAnc < -9), ]
a2.lowZ %>%
  tidyr::gather(., "pop", "anc", 1:14) %>%
  ggplot(., aes(x = pop, y = anc)) +
  geom_boxplot()
a2[which(a2.ZAnc < 0 & a2.ZAnc > -3.5), ] %>%
  tidyr::gather(., "pop", "anc", 1:14) %>%
  ggplot(., aes(x = pop, y = anc)) +
  geom_boxplot()
a2 %>%
  mutate(category_ZAnc = ifelse(a2.ZAnc < -3.5, "very low", 
                                ifelse(a2.ZAnc < 3.5 & a2.ZAnc >= -3.5, "middle", 
                                       ifelse(a2.ZAnc >=3.5 & a2.ZAnc < 10, "high", 
                                              "very  high")))) %>%
  tidyr::gather(., "pop", "anc", 1:14) %>%
  ggplot(., aes(x = pop, y = anc, color = category_ZAnc)) +
  geom_boxplot() +
  ggtitle("pop ancestries for very low to very high ZAnc scores")
# looks like in general pops have low ancestry freq (good) for very low ZAnc scores,
# BUT any population could still be a high outlier for the very lowest ZAnc scores
ggsave("plots/pop_anc_by_ZAnc_score.png", 
       width = 16, height = 8, units = "in",
       device = "png")

# ok can I replicated getting a ZAnc score < -3.5, or 0 mex anc for all pops?
# what do I get if I just set one population to freq 1 and all others to 0?
apply(diag(14), 2, function(x) ZAnc(ancFreq = x, invL = calcInvL(a2.cov), alpha = a2.alpha))
# yes I can get extreme negative values in this case..lower than all zeros
zAnc010 <- do.call(cbind,
                   lapply(1:14, function(i) 
                     zAnc(ancFreq = diag(14)[, i], invL = calcInvL(a2.cov), alpha = a2.alpha)))
apply(zAnc010, 2, sum)
zAnc(ancFreq = rep(0,14), invL = calcInvL(a2.cov), alpha = a2.alpha)
# or 0.5?
0.5*diag(14)
apply(0.5*diag(14), 2, function(x) ZAnc(ancFreq = x, invL = calcInvL(a2.cov), alpha = a2.alpha))

# what is the theoretical max if all pops have freq = 1? ~29; we're not hitting it
ZAnc(ancFreq = rep(1, length(a2.alpha)), invL = calcInvL(a2.cov), alpha = a2.alpha)
# what about if you have all but one population at very high freq, can we exceed ZAnc >29? Also, yes.
matrix101<- matrix(1,14,14)
diag(matrix101)<-0
zAnc101 <- do.call(cbind,
                   lapply(1:14, function(i) 
                     zAnc(ancFreq = matrix101[, i], invL = calcInvL(a2.cov), alpha = a2.alpha)))
apply(zAnc101, 2, sum)

# can I plot the cholesky matrix to look at what's going on?
melt(t(chol(a2.cov))) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  ggtitle("cholesky L - maize")
ggsave("plots/maize_chol_L.png", 
       width = 8, height = 8, units = "in",
       device = "png")

# and the inverse cholesky?
melt(solve(t(chol(a2.cov)))) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  ggtitle("inverse cholesky - maize")
ggsave("plots/maize_inverse_chol_L.png", 
       width = 8, height = 8, units = "in",
       device = "png")

    
#solve(t(chol(K)))  
ZAnc(ancFreq = rep(0, length(a2.alpha)), invL = solve(t(chol(a2.K))), alpha = a2.alpha)
ZAnc(ancFreq = rep(0, length(a2.alpha)), invL = solve((chol(a2.K))), alpha = a2.alpha)
a2.BAD.ZAnc <- apply(a2, 1, function(x) 
  ZAnc(ancFreq = x, invL = solve(chol(a2.K)), alpha = a2.alpha))
summary(a2.BAD.ZAnc)
plot(a2.BAD.ZAnc, apply(a2, 1, mean))

hist(apply(a2, 1, mean)) # mean across pops
pairs(a2[500:1000,7:14])
apply(a2, 2, var)

a2.cor <- round(cor(a2), 3)
# plot population-by-population 
# ancestry covariance matrix
melt(a2.cov) %>%
  #mutate(pop1 = a2.pops$LOCALITY) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/variance_covariance_matrix_maize.png", 
       width = 8, height = 8, units = "in",
       device = "png")

# plot correlation
melt(a2.cor) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
min(a2.cor)
# def. none are negative -- doesn't correct pop size
# in sample var/cor, but ok for now

# what variance would I get if it were just binomial draws?
# npq
ggplot(a2.pops, aes(x = var_npq, y = var_observed)) +
  geom_point(aes(color = LOCALITY)) +
  geom_abline() +
  xlim(0, .035) +
  ylim(0, .035)
ggsave("plots/variance_within_pop_vs_binomial.png", 
       width = 8, height = 8, units = "in",
       device = "png")

# what happens when it's centered?
# should make no difference
alpha <- apply(a2, 2, mean) # mean per pop
a2.c <- t(apply(a2, 1, function(row) row - alpha)) # centered
# plot centered covariances:
round(cov(a2.c), 3) %>%
  melt(.) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
round(cov(a2.c), 3) %>%
  min(.)
round(cor(a2), 3) %>%
  melt(.) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
round(cor(a2), 3) %>%
  min(.)

# can I simulate a distribution of ancestries 
# under the observed covariance matrix?
# and just under a case with variance but
# no covariance in ancestry?

# and how about just straight binomial draws?
hist(rbinom(n = 1000, size = 8, prob = .2))
# simulate under binomial, 10,000 draws
a2.binom <- t(sapply(1:nrow(a2.pops), function(x) rbinom(n=10000, size = (a2.pops$count*2)[x], prob = a2.pops$alpha[x])/(a2.pops$count*2)[x]))
a2.pops$var_sim_binom <- apply(a2.binom, 1, var)
ggplot(a2.pops, aes(x = var_npq, y = var_sim_binom)) +
  geom_point(aes(color = LOCALITY)) +
  geom_abline() # variance is what we'd expect, not correction needed

# under simple binomial, 
# how unlikely is it to get any one pop > 0.8?
# all pops under 0.01? -- pretty unlikely, but one issue
# is there's no measurement error here 
# so that really freq=0 w/ small N per pop
table(apply(a2.binom, 2, function(i) sum(i<.01)))
# probability for individual pops of having
# no mexicana ancestry at a locus is moderate:
summary(sapply(1:nrow(a2.pops), function(x) pbinom(0, size = (a2.pops$count*2)[x], prob = a2.pops$alpha[x])))
# but probability of all of them together having
# no mexicana ancestry is very low
# (this doesn't include variance from drift or hmm model inference error!)
# we have big multiple testing issues
# but still it looks like there's some enrichment
# because prob is ~ 10^-11 but about 0.001 
# proportion of sites have mex < .01 at all pops
prod(sapply(1:nrow(a2.pops), function(x) pbinom(0, size = (a2.pops$count*2)[x], prob = a2.pops$alpha[x])))

# are there cases where all populations appear to 
# have adaptive introgression?
table(apply(a2, 1, function(i) sum(i>.9)))
table(apply(a2, 1, function(i) sum(i>.8)))
table(apply(a2, 1, function(i) sum(i>.7)))
# threshold needs to be fairly low to get 
# any peaks including all pops
table(apply(a2, 1, function(i) sum(i>.55)))

# plot these results:
png(filename = "plots/mex_ancestry_shared_peaks.png",
    width = 15, height = 15, units = "in", 
    res = 300)
par(mfrow=c(2,2))
for (t in c(.9, .8, .7, .55)){
  hist(apply(a2[apply(a2, 1, function(x) sum(x>t)>=1), ], 1, function(i) sum(i>t)),
       xlab = "# pops outliers at peak",
       xlim = c(1,14),
       #ylim = c(0, 50000),
       main = paste0("peak > ", t, " mex anc"),
       col = "blue")
}
par(mfrow=c(1,1))
dev.off()

# probability under binomial of all > 50%:
# fairly low for individual pops (excluding idea of drift)
summary(sapply(1:nrow(a2.pops), function(x) pbinom((a2.pops$count)[x], size = (a2.pops$count*2)[x], prob = a2.pops$alpha[x], lower.tail = F)))
# very very low probability under binomial of all pops > 50% freq
prod(sapply(1:nrow(a2.pops), function(x) pbinom((a2.pops$count)[x], size = (a2.pops$count*2)[x], prob = a2.pops$alpha[x], lower.tail = F)))
# probability individual pop has 100% mex freq
# basically zero? check this
summary(1 - sapply(1:nrow(a2.pops), function(x) pbinom((a2.pops$count*2)[x], size = (a2.pops$count*2)[x], prob = a2.pops$alpha[x], lower.tail = T)))


# what about cases with just one population
# is under selection for mex. ancestry?
# let's say one pop has very high freq
# and all other pops have < 50%
table(apply(a2[apply(a2, 1, function(x) all(x >.9 | x < .5)), ], 
            1, function(i) sum(i>.9)))
table(apply(a2[apply(a2, 1, function(x) all(x >.9 | x < alpha)), ], 
            1, function(i) sum(i>.9)))
# maybe it's more clear to ask, in cases where at
# least some pops have very high ancestry freq
# how many pops are at or below their mean?
table(apply(a2[apply(a2, 1, function(x) sum(x >.9)>=1), ], 
            1, function(i) sum(i<=alpha)))
hist(apply(a2[apply(a2, 1, function(x) sum(x >.9)>=1), ], 
           1, function(i) sum(i<=.2)),
     main = "# pops w/ below 20% mex anc at peak",
     xlab = "peak > .9",
     freq = F)
# slightly lower threshold for a peak and there may be
# a few single pop selection events:
table(apply(a2[apply(a2, 1, function(x) sum(x >.8)>=1), ], 
            1, function(i) sum(i<=alpha)))
# save plots:
png(filename = "plots/mex_ancestry_non-shared_peaks.png",
    width = 15, height = 8, units = "in", 
    res = 300)
par(mfrow=c(1,3))
for (t in c(.9, .8, .6)){
  hist(apply(a2[apply(a2, 1, function(x) sum(x >t)>=1), ], 
             1, function(i) sum(i<=alpha)),
       xlab = "# pops < mean mex anc at peak",
       main = paste0("peak > ", t, " mex anc a pop"),
       freq = F)
}
par(mfrow=c(1,1))
dev.off()

# doing binomial comparison here doesn't make sense
# because basically no simulation results get that high for ind. pops in freq

table(apply(a2.binom[ , apply(a2.binom, 2, function(x) sum(x >.7)>=1)], 
            2, function(i) sum(i<=alpha)))


# are there some 'deserts' of mexicana ancestry?
table(apply(a2, 1, function(i) sum(i<.01)))
# yes, though by chance likelihood may
# not be that low (will need to check)
png(filename = "plots/mex_ancestry_deserts.png",
    width = 10, height = 8, units = "in", 
    res = 300)
par(mfrow=c(1,2))
hist(apply(a2.binom, 2, function(i) sum(i<.01)), 
     freq = F,
     xlim = c(0, 14),
     col = "yellow",
     main = "binomial simulation",
     xlab = "# pops w/ <.01 mex ancestry")
hist(apply(a2, 1, function(i) sum(i<.01)), 
     freq = F, 
     xlim = c(0, 14),
     col = "blue",
     main = "mexicana ancestry 'deserts'",
     xlab = "# pops w/ <.01 mex ancestry")
     #col = alpha("blue", .5),
     #add = T)
par(mfrow=c(1,1))
dev.off()

# very lenient is x < mean ancestry for that pop
table(apply(a2, 1, function(i) sum(i<alpha)))


# to do: simulate data under a multivariate normal
# using observed variances and covariances
mvn = mvrnorm(n=10000, # create some MVN data
                mu = alpha, 
                Sigma = a2.cov)
hist(mvn) #not normal - a mixture of correlated normals
summary(mvn)

# visualize how well the covariance and mean alpha were maintained
mvn.cov <- round(cov(mvn), 3)
# sanity check -- should preserve perfectly
# covariance matrix
melt(mvn.cov) %>%
  ggplot(data = ., aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
# and population means
plot(a2.pops$alpha ~ colMeans(mvn))
# simulate mvn without any covariances between pops
mvn.var.only = mvrnorm(n=10000, # create some MVN data
                       mu = alpha, 
                       Sigma = diag(a2.pops$var_observed))
hist(mvn.var.only)

# plot differences in mean of pop. mean ancestries
# for observed data and different model assumptions
png(file = "plots/data_mvn_binomial_mean_pop_anc.png",
    width = 15, height = 15, units = "in", res = 300)
par(mfrow = c(2,2))
hist(rowMeans(a2),
     xlim = c(0,1),
     main = "Obs. mean of pops",
     freq = F,
     col = "blue")
hist(rowMeans(mvn), 
     xlim = c(0, 1),
     main = "MVN mean of pops",
     col = "red",
     freq = F)
hist(rowMeans(mvn.var.only), 
     xlim = c(0, 1),
     main = "MVN (omit cov.) mean of pops",
     freq = F,
     col = "orange")
hist(colMeans(a2.binom),
     xlim = c(0,1),
     main = "Binom. mean of pops",
     freq = F,
     col = "yellow")
par(mfrow=c(1,1))
dev.off()

# what is my expected variance given the variance-covariance
# matrix observed?
# new variance for total = sum of variances plus sum of covariances
# then variance of the mean is divided by 14^2 
#var_exp_mean_pop_freq <- sum(a2.cov[upper.tri(a2.cov, diag = T)])/14
var_exp_mean_pop_freq <- sum(a2.cov)/14^2
mean_exp_mean_pop_freq <- mean(a2.pops$alpha)
# my variance is correct -- 
# the sum of correlated normals has variance equal to sum of variances & covariances
png(file = "plots/data_mvn_curve_mean_pop_anc.png",
    width = 15, height = 20, units = "in", res = 300)
par(mfrow = c(3,1))
hist(rowMeans(mvrnorm(n=1000000, # create a LOT of MVN data
                            mu = alpha, 
                            Sigma = a2.cov)), 
     #xlim = c(0, 1),
     #ylim = c(0, 5),
     xlim = c(-.1,1),
     ylim = c(0, 5.5),
     main = "MVN mean of pops",
     breaks = 100,
     freq = F)
abline(v = mean(rowMeans(mvn)), col = "black")
abline(v = mean_exp_mean_pop_freq, col = "blue")
curve(dnorm(x, 
            mean_exp_mean_pop_freq,
            sd = sqrt(var_exp_mean_pop_freq)),
      col = "red", add = T)
# the mvn is an ok approximation to the true data
hist(rowMeans(a2),
     xlim = c(-.1,1),
     ylim = c(0, 5.5),
     main = "Obs. mean of pops",
     freq = F,
     col = "blue")
curve(dnorm(x, # regular normal
             mean_exp_mean_pop_freq,
             sd = sqrt(var_exp_mean_pop_freq)),
      col = "red", add = T)
#curve(dcnorm(x, # censored normal
#            mean_exp_mean_pop_freq,
#            sd = sqrt(var_exp_mean_pop_freq),
#            left = 0),
#      col = "green", add = T)
hist(rcnorm(100000, # censored normal
            mean_exp_mean_pop_freq,
            sd = sqrt(var_exp_mean_pop_freq),
            left = 0),
     xlim = c(-.1,1),
     ylim = c(0, 5.5),
     col = "pink",
     main = "censored normal mean of pops",
     freq = F)
curve(dcnorm(x, # censored normal
             mean_exp_mean_pop_freq,
             sd = sqrt(var_exp_mean_pop_freq),
             left = 0),
      col = "black", add = T)
dev.off()

# this isn't quite the correct solution
# because the truncated normal doesn't preserve
# the data's variance and mean (oops .. maybe a correction exists)
a2.cnorm <- rcnorm(10000, # censored normal
                   mean_exp_mean_pop_freq,
                   sd = sqrt(var_exp_mean_pop_freq),
                   left = 0)

cov(mvn[,1], mvn[,2]) # empirical is close to expected .2
chol_mvn = t(chol(a2.cov)) # take the transpose to get a lower triangular matrix
chol_mvn %*% t(chol_mvn) == a2.cov #Sigma = LL'
  #inv_chol_mvn = solve(chol_mvn)
inv_chol_mvn = calcInvL(a2.cov)
z_mvn = apply(mvn, 1, function(i) inv_chol_mvn %*% (i - alpha))
cor(z_mvn[1,], z_mvn[2,]) # no correlation anymore, great
cor(z_mvn[3,], z_mvn[2,])
cor(z_mvn[3,], z_mvn[1,])
hist(z_mvn) # ~N(0,1)
hist(rnorm(50000*3, mean=0, sd=1), col="darkblue", add=TRUE)
hist(z_mvn[1,])# ~ N(0,1) 


# TO DO: I could look at totals across all pops
# rather than averaging within pop and then across pops
a <- a2[1,]
a2.alt <- apply(a2, 1, function(a) 
  t(calcInvL(cov(a2)) %*% (a - a2.alpha)) %*% calcInvL(cov(a2))%*%rep(1,14))


t(calcInvL(cov(a2)) %*% (rep(1,14) - a2.alpha)) %*% calcInvL(cov(a2))%*%rep(1,14)
t(calcInvL(cov(a2)) %*% (rep(0,14) - a2.alpha)) %*% calcInvL(cov(a2))%*%rep(1,14)
t(calcInvL(cov(a2)) %*% (a2.alpha - a2.alpha)) %*% calcInvL(cov(a2))%*%rep(1,14)

a2.invL <- calcInvL(cov(a2))
a2.alt <- t(apply(a2, 1, function(a) 
  t(a2.invL %*% (a - a2.alpha)) %*% (a2.invL%*%rep(1,14))))
plot(apply(a2, 1, mean), a2.alt, col = "blue")
abline(h = t(calcInvL(cov(a2)) %*% (rep(1,14) - a2.alpha)) %*% calcInvL(cov(a2))%*%rep(1,14)
, col = "green")
abline(h = t(calcInvL(cov(a2)) %*% (rep(0,14) - a2.alpha)) %*% calcInvL(cov(a2))%*%rep(1,14)
       , col = "orange")

t(apply(diag(14), 1, function(x) 
  t(a2.invL %*% (x - a2.alpha)) %*% (a2.invL%*%rep(1,14))))
t(apply(diag(14), 1, function(x) 
  t(a2.invL %*% (x - a2.alpha)) %*% (a2.invL%*%rep(-1,14))))


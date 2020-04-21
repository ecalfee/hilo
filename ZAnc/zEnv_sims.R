# simple simulations of environmental associations + MVN noise
library(dplyr)
library(tidyr)
library(ggplot2)
source("k_matrix.R")

# functions
# truncate to [0,1] range
truncate01 <- function(x){
  t <- x
  t[t < 0] <- 0
  t[t > 1] <- 1
  return(t)
}


# read in data (or could pick arbitrary alpha's and K matrices)
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

maize_pops <- unique(meta.pops$pop[meta.pops$zea == "maize"])


K <- read.table("results/zAnc_K_matrix_maize.txt", sep = "\t", header = T) %>%
  as.matrix(.)
alpha <- read.table("results/zAnc_alpha_maize.txt", sep = "\t", header = T)$x
names(alpha) <- colnames(K)
n_sim = 1000
set.seed(101)
# mvn (neutral)
mvn <- mvrnorm(n = n_sim, 
                         mu = alpha, 
                         Sigma = K, 
                         tol = 1e-6, 
                         empirical = FALSE, 
                         EISPACK = FALSE)
# add additional slope with elevation (units = km)
elev = sapply(names(alpha), function(p)
              meta.pops$ELEVATION[meta.pops$pop == p]/1000)
elev_c = elev - mean(elev)
slopes = do.call(c, lapply(c(0, 0.05, 0.5, 5, -0.05, -0.5, -5), rep, 3))
intercepts = rep(c(0, 0.3, -0.3), 7)
sel_names = 1:length(slopes)
lapply(1:length(slopes), function(i) summary(t(alpha + intercepts[i] + elev_c*slopes[i])))
sel <- lapply(1:length(slopes), function(i)
  t(apply(mvn, 1, function(x) x + intercepts[i] + elev_c*slopes[i])))

# truncate 0-1 range
mvn_trunc <- truncate01(mvn)
sel_trunc <- lapply(sel, truncate01)

# what I really want to do is not a logistic, but add slope, then truncation at 0, 1
# but logistic might be a good alternative. apply logistic before or after noise?
X = cbind(intercept = 1, elev_c) %>%
  as.matrix(.)
invK = solve(K)
#y1 = mvn[1,] # first example
#invK %*% K %>%
#  round(., 5)
# generalized least squares model for population ancestries centered at pop mean ancestry:
# (y - alpha)
mvn_betas = t(apply(mvn, 1, function(y) 
  solve(t(X) %*% invK %*% X) %*% t(X) %*% invK %*% (y - alpha)))
sel_betas = lapply(sel, function(s)
  t(apply(s, 1, function(y) 
    solve(t(X) %*% invK %*% X) %*% t(X) %*% invK %*% (y - alpha))))
mvn_trunc_betas = t(apply(mvn_trunc, 1, function(y) 
  solve(t(X) %*% invK %*% X) %*% t(X) %*% invK %*% (y - alpha)))
sel_trunc_betas = lapply(sel_trunc, function(s)
  t(apply(s, 1, function(y) 
    solve(t(X) %*% invK %*% X) %*% t(X) %*% invK %*% (y - alpha))))
# plot -- how well did we recover true values without truncation?
# with truncation?
# how did alphas fit in? Or do I need to subtract them off first..
# what is my power in the diff selection schemes to be in empirical top 1% neutral?

mvn_betas %>%
  as.data.frame(.) %>%
  data.table::setnames(c("b0", "bElev")) %>%
  pivot_longer(., colnames(.), names_to = "b", values_to = "estimate") %>%
  ggplot(., aes(x = estimate, fill = b)) +
  geom_histogram() +
  facet_wrap(~b, scales = "free_x") +
  ggtitle("MVN null model - no truncation")
mvn_betas %>%
  as.data.frame(.) %>%
  data.table::setnames(c("b0", "bElev")) %>%
  #pivot_longer(., colnames(.), names_to = "b", values_to = "estimate") %>%
  ggplot(., aes(x = b0, y = bElev)) +
  geom_point() +
  ggtitle("MVN null model - no truncation")

# with truncation:
mvn_trunc_betas %>%
  as.data.frame(.) %>%
  data.table::setnames(c("b0", "bElev")) %>%
  pivot_longer(., colnames(.), names_to = "b", values_to = "estimate") %>%
  ggplot(., aes(x = estimate, fill = b)) +
  geom_histogram() +
  facet_wrap(~b, scales = "free_x") +
  ggtitle("MVN null model - [0,1]")
mvn_trunc_betas %>%
  as.data.frame(.) %>%
  data.table::setnames(c("b0", "bElev")) %>%
  ggplot(., aes(x = b0, y = bElev)) +
  geom_point() +
  ggtitle("MVN null model - [0,1]")

# now show distribution of estimates for selection scenarios
sel_betas_pretty <- do.call(rbind, 
        lapply(1:length(sel_betas), function(i)
          as.data.frame(sel_betas[[i]]) %>%
            data.table::setnames(c("b0", "bElev")) %>%
            pivot_longer(., colnames(.), names_to = "b", values_to = "estimate") %>%
            mutate(name = sel_names[i],
                   intercept = intercepts[i],
                   slope = slopes[i])))

# add lines for truth. also just add in the null model at the top
sel_betas_pretty %>% # top row is the neutral model
  ggplot(., aes(x = estimate, fill = b)) +
  geom_histogram(bins = 100) +
  geom_vline(data = data.frame(name = sel_names,
                               b0 = intercepts,
                               bElev = slopes) %>%
               pivot_longer(., c("b0", "bElev"), names_to = "b", values_to = "estimate"),
               aes(xintercept = estimate)) +
  facet_grid(name ~ b, scales = "free_x")

# power to detect a selected locus at 5% FDR:
#table(mvn == sel[[1]]) # all True
# start with extreme scenario 10:
filter(sel_betas_pretty, b == "bElev") %>%
  group_by(name) %>%
  summarise(max = max(estimate))
head(sel_betas_pretty$estimate)
vals_bElev = seq(from = 0, to = 6, by = .01)
over_null = sapply(vals_bElev, function(x) sum(sel_betas[[1]][ , 2] > x))
over_10 = sapply(vals_bElev, function(x) sum(sel_betas[[10]][ , 2] > x))
# assume 1% of loci are selected, what's your false-discovery-rate cutoff?
perc_sel = .01
FDR_10 = over_null/(over_null + perc_sel * over_10)
over_10[min(which(FDR_10 < .05))]/n_sim # 100% power
pwr <- function(null_b, sel_b, perc_sel, v, fdr = 0.05){
  # input: null model results, selection model results, percent of loci selected
  # v vector of test values to evaluate fdr at, 
  # fdr threshold to evaluate power to find the selected site
  over_null = sapply(v, function(x) sum(null_b > x))/length(null_b)
  over_sel = sapply(v, function(x) sum(sel_b > x))/length(sel_b)
  fdrs = over_null/(over_null + perc_sel * over_sel)
  p = over_sel[min(which(fdrs < fdr))]
  return(p)
}
pwr(null_b = sel_betas[[1]][ , 2], sel_b = sel_betas[[10]][ , 2],
    perc_sel = 0.01, v = vals_bElev, fdr = 0.05)
sapply(sel_betas, function(s) pwr(null_b = sel_betas[[1]][ , 2],
                                  sel_b = s[ , 2],
                                  perc_sel = 0.01,
                                  v = vals_bElev,
                                  fdr = 0.05))
sapply(sel_trunc_betas, function(s) pwr(null_b = sel_trunc_betas[[1]][ , 2],
                                  sel_b = s[ , 2],
                                  perc_sel = 0.01,
                                  v = vals_bElev,
                                  fdr = 0.05))

# now look at the real data. what are the b0 and bElev estimates genomewide and in inv4m?

# now with truncation:
sel_trunc_betas_pretty <- do.call(rbind, 
                            lapply(1:length(sel_trunc_betas), function(i)
                              as.data.frame(sel_trunc_betas[[i]]) %>%
                                data.table::setnames(c("b0", "bElev")) %>%
                                pivot_longer(., colnames(.), names_to = "b", values_to = "estimate") %>%
                                mutate(name = sel_names[i],
                                       intercept = intercepts[i],
                                       slope = slopes[i])))

# add lines for truth. also just add in the null model at the top
sel_trunc_betas_pretty %>% # top row is the neutral model
  ggplot(., aes(x = estimate, fill = b)) +
  geom_histogram(bins = 100) +
  geom_vline(data = data.frame(name = sel_names,
                               b0 = intercepts,
                               bElev = slopes) %>%
               pivot_longer(., c("b0", "bElev"), names_to = "b", values_to = "estimate"),
             aes(xintercept = estimate)) +
  facet_grid(name ~ b, scales = "free_x")




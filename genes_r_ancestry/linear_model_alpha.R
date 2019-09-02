# in this script I fit a linear model for the admixture proportion alpha
# using estimated allele frequencies for allopatric and sympatric populations

library(dplyr)
library(rethinking)
# use multiple cores for model fits
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

pops = c("allo.maize", "allo.mexicana", "symp.maize", "symp.mexicana", "parv")
# paths to input allele frequencies
path_allele_freqs = "../within_ancestry_diversity/results/allele_freq/pass2_alloMAIZE"
# region
region=14
# sites file
sites_file = paste0("../variant_sites/results/pass2_alloMAIZE/region_", region, ".var.sites")

sites0 <- read.table(sites_file, 
                     stringsAsFactors = F, header = F) %>%
  data.table::setnames(c("scaffold", "pos", "major", "minor"))

allele_freq <- lapply(pops, function(POP)
                      read.table(paste0(path_allele_freqs, "/", POP, "/region_", region, ".mafs.gz"),
                          stringsAsFactors = F, header = T) %>%
  left_join(sites0, ., by = c("scaffold"="chromo", "pos"="position",
                              "major", "minor")) %>%
    dplyr::select(scaffold, pos, phat) %>%
    data.table::setnames(c("scaffold", "pos", POP)))
f <- bind_cols(sites0,
               do.call(cbind, lapply(pops, function(POP)
  read.table(paste0(path_allele_freqs, "/", POP, "/region_", region, ".mafs.gz"),
             stringsAsFactors = F, header = T) %>%
    left_join(sites0, ., by = c("scaffold"="chromo", "pos"="position",
                                "major", "minor")) %>%
    dplyr::select(phat) %>%
    mutate(phat = ifelse(phat > 1, 1, phat)))) %>% # rounding error in angsd? estimated allele freqs very slightly > 1
  data.table::setnames(pops))
with(f, lm(symp.maize ~ allo.maize + allo.mexicana))
with(f, lm(symp.mexicana ~ allo.maize + allo.mexicana))
with(f, lm(parv ~ allo.maize + allo.mexicana))

m_maize <- map( # quadratic approximation of the posterior MAP
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
  data = f[complete.cases(f), ],
  start = list(theta = 2, alpha = 0.5))
precis(m_maize)
pairs(m_maize)
# stan version
m2_maize <- map2stan( # quadratic approximation of the posterior MAP
  alist(
    #symp.maize ~ dnorm(p, theta),
    symp.maize ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    p <- alpha*allo.mexicana + (1 - alpha)*allo.maize,
    alpha ~ dunif(0, 1),
    theta ~ dunif(0, 10)
  ),
  # needs no NAs and needs to be a SNP segregating between allopatric maize and allopatric mexicana
  data = f[complete.cases(f) & !(f$allo.maize == 0 & f$allo.mexicana == 0), ],
  start = list(theta = 2, alpha = 0.5),
  chains = 2, iter = 1000, warmup = 500)
precis(m2_maize)
pairs(m2_maize)



m_mexicana <- map( # quadratic approximation of the posterior MAP
  alist(
    symp.mexicana ~ dnorm(p, theta),
    #symp.maize ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    p <- alpha*allo.mexicana + (1 - alpha)*allo.maize,
    alpha ~ dunif(0, 1),
    theta ~ dunif(0, 10)
  ),
  data = f[complete.cases(f), ],
  start = list(theta = 2, alpha = 0.5))
precis(m_mexicana)
pairs(m_mexicana)

m_parv <- map( # quadratic approximation of the posterior MAP
  alist(
    parv ~ dnorm(p, theta),
    #symp.maize ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    p <- alpha*allo.mexicana + (1 - alpha)*allo.maize,
    alpha ~ dunif(0, 1),
    theta ~ dunif(0, 10)
  ),
  data = f[complete.cases(f), ],
  start = list(theta = 2, alpha = 0.5))
precis(m_parv)
pairs(m_parv)


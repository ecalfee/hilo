# this script has helper functions to
# flatten or 'melt' covariance matrices K
# so they can be plotted
# population by population covariance matrix K:

# (I need haps_per_pop so I can do a finite sample size correction for
# the variance within pop to fairly compare it with
# covariances between populations)
kmelt <- function(K, haps_per_pop, correlation = F, regions = NULL){ 
  # haps per pop is to do the finite sample size correction
  # for the variance (on the diagonal) so it can be compared
  # to covariances between populations
  
  # want matrix 1's everywhere else but (N-1)/N on diagonals--
  K_corrected = K*(matrix(1, dim(K)[1], dim(K)[2]) + diag((haps_per_pop-1)/haps_per_pop))
  # alternatively look at correlations, not covariances
  if (correlation) K_corrected = cov2cor(K_corrected)
  k_melted = reshape2::melt(K_corrected, 
                            value.name = "covariance")[upper.tri(K_corrected, diag = T),] %>%
    separate(., Var1, c("zea_SHORT_1", "LOC_SHORT_1")) %>%
    separate(., Var2, c("zea_SHORT_2", "LOC_SHORT_2")) %>%
    mutate(., same_location = LOC_SHORT_1 == LOC_SHORT_2) %>%
    mutate(., mex_vs_maize = zea_SHORT_1 != zea_SHORT_2) %>%
    mutate(., have_pairs = LOC_SHORT_1 %in% pop_pairs_SHORT & LOC_SHORT_2 %in% pop_pairs_SHORT)
  if (! is.null(regions)){ # label these results
    k_melted$regions <- regions
  }
  return(k_melted)
}

# individual by individual covariance matrix K
kmelt_ind <- function(K, regions = NULL, correlation = F){ 
  # with individuals, the finite sample size correction
  # for the variance (on the diagonal) is N = 2, 
  # so it can be compared to individual covariances 
  # between populations
  
  # want matrix 1's everywhere else but N-1/N on diagonals = 1/2 (b/c diploid) --
  K_corrected = K*(matrix(1, dim(K)[1], dim(K)[2]) + diag(rep(1/2, dim(K)[1])))
  if (correlation) K_corrected = cov2cor(K_corrected)
  k_melted = reshape2::melt(K_corrected, 
                            value.name = "covariance")[upper.tri(K_corrected, diag = T),] %>%
    mutate(., ID_1 = as.character(Var1)) %>%
    mutate(., ID_2 = as.character(Var2)) %>%
    left_join(., meta, by = c("ID_1"="ID")) %>%
    left_join(., meta, by = c("ID_2"="ID"), suffix = c("_1", "_2")) %>%
    mutate(., same_location = LOCALITY_1 == LOCALITY_2) %>%
    mutate(., mex_vs_maize = zea_1 != zea_2) %>%
    mutate(., same_ind = ID_1 == ID_2) %>%
    mutate(., have_pairs = (LOCALITY_1 %in% pop_pairs && LOCALITY_2 %in% pop_pairs))
  if (! is.null(regions)){ # label these results
    k_melted$regions <- regions
  }
  return(k_melted)
}
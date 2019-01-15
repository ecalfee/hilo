library(dplyr)
library(tidyr)
library(ggplot2)
# Calculating Ast, analogous to Fst
# Hudson 1995 estimator of Fst using sample allele frequencies
# p1 and p2 with sample sizes n1 and n2 (eq. 10 in Bhatia 2013)
# Here I'm assuming that the variance of the allele frequency estimator
# is binomial, where n is number of haplotypes sampled (see derivation in Bhatia 2013 supplement)
hudson_95_fst <- function(p1, p2, n1, n2){
  het_within <- (p1-p2)^2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
  het_total <- p1*(1-p2) + p2*(1-p1)
  fst <- het_within/het_total
  return(fst)
}
# function to calculate Ast for each pair of populations (return flattened into a dataframe)

# calculate for all, low_r and high_r recombination rate bins

# plot as a covariance matrix

# plot as boxplots:
# (1) comparing within and between locations 
# for mexicana-maize pairs

# (2) maize pairs compared to mexicana pairs 
# compared to maize-mex pairs

# (3) high vs. low vs. all recombination rates, 
# facet wrap by w/ and w/out inversion and maize vs. mex
# (ignore maize-mex pairs)
 
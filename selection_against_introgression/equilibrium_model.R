library(data.table)
library(dplyr)

# script to estimate initial admixture proportion and mean
# selection against introgresed alleles: s, mu, alpha0

# TO DO: load in some data
# handle special cases (no SNPs in focal region to calculate f4, 
# no exons, smaller window due to edge of chromosome), very large exons

# maximum distance (in cM) to look for 'nearby' exons
# around the center of a focal region
max_dist_cM <- 1 # so makes 2cM windows

# grid s-values for mean selection against del. introgressed variants
# and mu-values for mean per bp occurence of deleterious variants
grid_s <- c(10^-6, 10^-5, 10^-4, 10^-3)
grid_u <- c(10^-7, 10^-6, 10^-5, 10^-4)

# observed f4 statistics, in 1kb windows (or .01cM -- does this get too small?)
f4_admix <- c(.2, .3, .4, .1, .3)
f4_nonadmix <- c(.4, .5, .4, .2, .5)

# exon chr, center positions (cM_pos), and length (bp)
coding <- data.frame(
  chr = rep(1, 10), # in order
  cM_pos = cumsum(runif(10, 0, 1)), # in order
  coding_length = c(100, 200, 300, 200, 250, 80, 300, 500, 250, 10)
)

# focal region positions (cM), and chromosome
pos <- cumsum(runif(10, 0, 1))[c(1,3,5,7,9)] # just make 5 random window centers
chr <- rep(1, 5)

# make a list of dataframes,
# one dataframe per focal position,
# saving the nearby exons (done only once)
list_nearby_exons <- lapply(1:length(pos), function(i)
                           filter(coding, chr == chr[i]) %>%
                             mutate(., dist_cM = abs(cM_pos - pos[i])) %>%
                             filter(., dist_cM < max_dist_cM) %>%
                             dplyr::select(., c("dist_cM", "coding_length")) %>% # delete extra columns to save space
                             dplyr::arrange(., dist_cM))

# calculate expected equilibrium fraction of ancestry surviving selection
# for a given s, u, and focal region, using two vectors holding information 
# about each coding feature within a larger window around focal region:
# dist_M = (Morgans) distance from focal region to center of each coding feature
# coding_length = length in bp for each coding feature
# note: coding features should be ordered by small -> large distance from focal region
calc_g <- function(s, u, dist_cM, coding_length){
  dist_M <- dist_cM/100 # convert cM to Morgans
  m <- length(coding_length) # number of coding features in larger window
  prod(1 - u*coding_length) +   # 1 times the probability of no mutations in any features
    sum(dist_M/(dist_M + s) * # summing over each coding feature j, r/(r+s) times
          u*coding_length * # probability of a mutation in coding feature j
          data.table::shift(x = cumprod(1 - u*coding_length), n = 1, fill = 1)) # and none closer
}

# calculate alpha0; analytical solution minimizing SDD(s,u,alpha0)
calc_alpha0 <- function(f4_admix, f4_nonadmix, g) {
  sum(f4_admix*g*f4_nonadmix)/sum((g*f4_nonadmix)^2)
}

# calculate SDD(s,u,alpha)
calc_SDD <- function(s, u, alpha0){
  sum((f4_admix - g*alpha0*f4_nonadmix)^2) # sum of squared diff measured f4 and model equilibrium prediction for f4
}

# grid search to minimize SSD(theta)
# make empty grid to fill in results
SDD_grid <- matrix(NA, length(grid_s), length(grid_u)) 
alpha0_grid <- matrix(NA, length(grid_s), length(grid_u))

# calculate g(s, u) for all focal regions
# g <- sappy(list_nearby_exons, function(x) 
#  calc_g(s = s, u = u, dist_M = x$))
for (i in 1:length(grid_s)){
  for (j in 1:length(grid_u)){
    s <- grid_s[i]
    u <- grid_u[j]
    g <- sapply(list_nearby_exons, function(x) # most of the computation happens here
                 calc_g(s = s, u = u, 
                        dist_cM = x$dist_cM, 
                        coding_length = x$coding_length))
    # calculate optimum alpha0 that minimizes SDD
    alpha0 <- calc_alpha0(f4_admix = f4_admix, 
                          f4_nonadmix = f4_nonadmix, 
                          g = g)
    # save optimum alpha0
    alpha0_grid[i, j] <- alpha0 
    # calculate SDD(s,u,alpha0) & save result
    SDD_grid[i, j] <- calc_SDD(s = s, u = u, alpha0 = alpha0)
  }
}



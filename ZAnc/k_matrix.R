#Kinship matrix using empirical pairwise ancestry covariance 
#(using individual genome-wide ancestry proportion alpha[i] as the 'mean' from which ancestry deviates)
#k[ij] = mean over loci l of (x[il] - alpha[i])*(x[jl] - alpha[j])

suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(poibin))
# for poisson binomial distribution

make_dip <- function(ms_file){
  #d <- readLines(ms_file, n=200) #max lines to read is set to 100 to not overwhelm memory - this should be a short sample
  #d <- d[-c(1, 2, 3, length(d))] #first 3 lines are header and last line is blank; remove
  #sample_size <- length(d)/2
  #n_markers <- nchar(d[1]) #need to stringsplit by character
  #file lists genotypes as continuous string of 111000111. The next line of code splits this into matrix entries of raw genotypes, 2 per diploid ind.
  #hap <- matrix(unlist(lapply(strsplit(d,""), as.numeric)), ncol = n_markers, byrow = T)
  #create matrix of haplotypes
  hap <- make_hap(ms_file)
  #create new matrix with half the rows for diploid (not haploid) genotypes
  sample_size <- dim(hap)[1]/2
  n_markers <- dim(hap)[2]
  dip <- matrix(nrow = sample_size, ncol = n_markers) 
  for (k in 1:sample_size) {
    dip[k,] <- hap[(2*k-1),] + hap[(2*k),] #each diploid individual has originally 2 rows that I condense into one by adding genotypes
  }
  return(dip)
}

make_hap <- function(ms_file, keep = TRUE){ # keep is an optional boolean vector to say which samples/rows to keep
  d <- readLines(ms_file)
  d <- d[-c(1, 2, 3, length(d))] #first 3 lines are header and last line is blank; remove
  d <- d[keep] #subsample as specified by keep, default = T (keep all)
  n_markers <- nchar(d[1]) #need to stringsplit by character
  #file lists genotypes as continuous string of 111000111. The next line of code splits this into matrix entries of raw genotypes, 2 per diploid ind.
  hap <- matrix(unlist(lapply(strsplit(d,""), as.numeric)), ncol = n_markers, byrow = T)
  return(hap)
}



plot_anc = function(anc, col0 = "#99FF66", col1 = "#009900", ...){ #plots haploid ancestry with dark and light green
  cols <- c(
    '0' = col0,
    '1' = col1
  )
  image(1:dim(anc)[2], 1:dim(anc)[1], t(anc), col=cols, main="ancestry", ...)
}

calcK = function(ancFreqMatrix, alpha){ # for all markers and # pops = length(alpha)
  ((ancFreqMatrix - alpha) %*% t(ancFreqMatrix - alpha))/ncol(ancFreqMatrix) #matrix multiply
  #above matrix operation is equivalent to the original code below
  #K <- matrix(nrow = sample_size, ncol = sample_size)
  #for (i in 1:sample_size) {
  #for (j in 1:sample_size) {
  #K[i,j] <- mean((anc[i,] - alpha[i])*(anc[j,] - alpha[j]))
  #}
  #}
}

calc_loss_anc_het = function(ancFreqMatrix, alpha){ # for all markers and # pops = length(alpha)
  #het_obs = ((ancFreqMatrix %*% t(1 - ancFreqMatrix)) + ((1 - ancFreqMatrix)*t(ancFreqMatrix)))/ncol(ancFreqMatrix) #matrix multiply
  #het_exp = t(alpha)*(1-alpha) + t(1-alpha)*alpha
  #above matrix operation is equivalent to the original code below
  sample_size = ncol(ancFreqMatrix)
  het_obs = matrix(nrow = sample_size, ncol = sample_size)
  het_exp = matrix(nrow = sample_size, ncol = sample_size)
  for (i in 1:sample_size) {
  for (j in 1:sample_size) {
  het_obs[i, j] <- mean(ancFreqMatrix[ , i]*(1 - ancFreqMatrix[ , j]) + (ancFreqMatrix[ , j]*(1 - ancFreqMatrix[ , i])))
  het_exp[i, j] <- alpha[i]*(1 - alpha[j]) + alpha[j]*(1 - alpha[i])                 
  }
  }
  return(list(obs = het_obs, exp = het_exp))
}
# test:
calc_loss_anc_het(ancFreqMatrix = rbind(c(.3,.4),c(.7,.8),c(.5,.6)), alpha = c(.5, .6))

calcInvL = function(K){
  solve(t(chol(K)))  
  #calculate the inverse Cholesky. t(chol(K)) is a lower left triangular matrix such that LL' = K
}

# wrapper functions to calculating alpha, K and inverse Cholesky of K
make_K_matrix <- function(anc, ploidy = 1){ # anc is assumed haploid by default
  # but can be diploid or a population-wide ploidy (like 10, 40-ploid)
  # passing an ancestry frequency, instead of counts, and ploidy = 1 is the same as providing counts and the correct ploidy
  sample_size <- dim(anc)[1]
  n_markers <- dim(anc)[2]
  
  # ensure K will be invertible:
  # take out genotypes with no variation (and therefore no covariation)
  no_vary <- apply(anc, 1, function(i) Reduce( "&", i == mean(i)))
  # and identical genotypes
  copy2 <- duplicated(anc)
  # rows to keep
  keep_id <- !(copy2 | no_vary)
  # apply to ancestry
  ancFrq <- anc[keep_id, ]/ploidy # divide by ploidy to get an ancestry freq from a count
  
  #individual genomewide ancestry proportions or global ancestry proportions
  alpha = apply(ancFrq, 1, mean) # to calc deviation of local ancestry from genomewide average, alpha
  
  #K covariance in ancestry deviation matrix
  K = calcK(ancFreqMatrix, alpha) 
  
#return(K)  
  invC = calcInvL(K) #calculate the inverse Cholesky. t(chol(K)) is a lower left triangular matrix such that LL' = K
  
  return(list(K=K, alpha=alpha, invC=invC, keep_id=keep_id)) # keep_id stores which rows to keep
  }

# function for calculating statistic:
make_K_calcs = function(pop_anc){
  # calculate mean population frequency across all snps
  pop_alpha = apply(pop_anc, 1, mean)
  # calculate K matrix for populations
  pop_K = calcK(ancFreqMatrix = pop_anc, alpha = pop_alpha)
  pop_InvL = calcInvL(pop_K)
  return(list(alpha = pop_alpha, K = pop_K, InvL = pop_InvL))
}


plot_corr_matrix <- function(ms_file, hap = T){# haploid option instead of diploid option (haploid by default)
  if (hap) corrplot(cov2cor(make_K_matrix(make_hap(ms_file))$K), method="shade")#, main = "haploid ancestry corr")
  else corrplot(cov2cor(make_K_matrix(make_dip(ms_file))$K), method="shade")#, main = "diploid ancestry corr")
}

# function to make a K matrix for either one locus or a pair of loci
make_K_d <- function(alpha, l1, l2){ # l1 and l2 are ancestries at 2 loci (either diploid or haploid ancestry) for all sampled ind.
  # if l1 = l2 it's a variance calculation within individuals and covariance for the same locus across individuals
  # if l1 != l2 it's a covariance or LD between loci calculation within individuals and a cross-individual covariance at diff loci across ind.
  # l1 and l2 can be vectors for one locus (length = # sampled inviduals) or matrices of n pairs of loci l1 and l2 recomb. distance d apart
  # in the special case where d = 0 and n = all sampled loci, this gives the K matrix above from make_K_matrix()

  if (is.vector(l1)) { # single locus
    n_markers = 1
  }
  
  else { # many pairs of loci
    n_markers = dim(l1)[2]
  }
  
  #K covariance in ancestry deviation matrix
  K = ((l1 - alpha) %*% t(l2 - alpha))/n_markers #matrix multiply
  # s.t. K[i,j] = mean covariance in ancestry deviation at loci l1 from individual i and loci l2 from individual j, distance d apart

  return(K)
}


# function to make a zTz matrix for either one locus or a pair of loci
make_zTz_matrix <- function(alpha, l1, l2, invC){ # l1 and l2 are ancestries at 2 loci (either diploid or haploid ancestry) for all sampled ind.
  # if l1 = l2 it's a variance calculation within individuals and covariance for the same locus across individuals
  # if l1 != l2 it's a covariance or LD between loci calculation within individuals and a cross-individual covariance at diff loci across ind.
  # l1 and l2 can be vectors for one locus (length = # sampled inviduals) or matrices of n pairs of loci l1 and l2 recomb. distance d apart
  # rotate by the cholesky constructed from a large number of neutral markers.
  
  if (is.vector(l1)) { # single locus
    n_markers = 1
  }
  
  else { # many pairs of loci
    n_markers = dim(l1)[2]
  }
  
  #K covariance in ancestry deviation matrix
  zTz_matrix = ((invC %*% (l1 - alpha)) %*% t(invC %*% (l2 - alpha)))/n_markers #matrix multiply
  # s.t. zTz[i,j] = mean standardized (w/ cholesky) covariance in ancestry deviation 
  # at loci l1 from individual i and loci l2 from individual j, distance d apart
  
  return(zTz_matrix)
}



if (FALSE){#silencing large block of code to build from later
  plot_anc(make_hap("10chromB/anc_haplotypes/ms_sim_output/sample40/panmixia/gen50/chr1_anc_1.txt")[, 5000:6000])
  plot_anc(make_hap("10chromB/anc_haplotypes/ms_sim_output/sample40/cont_mig_low/gen50/chr1_anc_18.txt"))
  plot_anc(make_hap("10chromB/anc_haplotypes/ms_sim_output/sample40/iso4/gen50/chr1_anc_18.txt"))
  plot_anc(make_hap("10chromB/anc_haplotypes/ms_sim_output/sample40/pairs_low/gen50/chr1_anc_8.txt"))
  plot_anc(make_hap("10chromB/anc_haplotypes/ms_sim_output/sample40/pairs_low/gen50/chr5_anc_8.txt"))
  plot_anc(make_hap("10chromB/anc_haplotypes/ms_sim_output/sample40/inv_low_pos/gen50/chr1_anc_8.txt"))
  plot_anc(make_hap("10chromB/anc_haplotypes/ms_sim_output/sample40/inv_low_neg/gen50/chr10_anc_8.txt"))  
  plot_anc(make_hap("10chromB/anc_haplotypes/ms_sim_output/sample40/inv_low_pos/gen50/chr1_anc_2.txt"))  
  plot_anc(make_hap("10chromB/anc_haplotypes/ms_sim_output/sample40/step4_low/gen50/chr1_anc_2.txt")) 
  #I'm not sure why I was plotting with make_dip diploid genotypes below. Clearly I should plot haplotypes
  plot_corr_matrix("10chromB/anc_haplotypes/ms_sim_output/sample40/panmixia/gen50/chr1_anc_1.txt")
  plot_corr_matrix("10chromB/anc_haplotypes/ms_sim_output/sample40/panmixia/gen50/chr1_anc_18.txt")
  plot_corr_matrix("10chromB/anc_haplotypes/ms_sim_output/sample40/panmixia/gen50/chr1_anc_21.txt")
  plot_corr_matrix("10chromB/anc_haplotypes/ms_sim_output/sample40/panmixia/gen50/chr7_anc_21.txt")
  plot_corr_matrix("10chromB/anc_haplotypes/ms_sim_output/sample40/cont_mig_low/gen50/chr7_anc_21.txt")
  plot_corr_matrix("10chromB/anc_haplotypes/ms_sim_output/sample40/pairs_low/gen50/chr1_anc_13.txt")
  plot_corr_matrix("10chromB/anc_haplotypes/ms_sim_output/sample40/iso4/gen50/chr1_anc_13.txt")
  plot_corr_matrix("10chromB/anc_haplotypes/ms_sim_output/sample40/iso4/gen50/chr5_anc_13.txt")
  plot_corr_matrix("10chromB/anc_haplotypes/ms_sim_output/sample40/inv_low_pos/gen50/chr5_anc_1.txt")
  plot_corr_matrix("10chromB/anc_haplotypes/ms_sim_output/sample40/inv_low_neg/gen50/chr1_anc_1.txt")
  plot_corr_matrix("10chromB/anc_haplotypes/ms_sim_output/sample40/inv_low_neg/gen50/chr10_anc_1.txt")
  
  plot_anc(make_dip("3chrom/anc_haplotypes/ms_sim_output/panmixia_anc/sample10_50gen_anc_1_chr1.txt"))
  plot_anc(make_dip("3chrom/anc_haplotypes/ms_sim_output/panmixia_anc/sample10_50gen_anc_1_chr3.txt"))
  plot_anc(make_dip("3chrom/anc_haplotypes/ms_sim_output/cont_mig_low_anc/sample10_50gen_anc_1_chr1.txt"))
  plot_anc(make_dip("3chrom/anc_haplotypes/ms_sim_output/panmixia_anc/sample10_50gen_anc_1_chr2.txt"))
  plot_anc(make_dip("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/step4_high_neutral/sample10_250gen_anc.txt"))
  plot_anc(make_dip("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/step4_high_1/sample10_250gen_anc.txt"))
  plot_anc(make_dip("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/cont_mig_low_1/sample10_250gen_anc.txt"))
  plot_anc(make_dip("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/inv_low_pos1/sample10_250gen_anc.txt"), xlab = "chromosome", ylab = "individuals") #good image to show for invasion with a selected allele
  plot_anc(make_dip("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/inv_low_neg1/sample10_250gen_anc.txt"), xlab = "chromosome", ylab = "individuals") #invasion with negatively selected alleles, also nice image
  plot_anc(make_dip("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/inv_high_neg1/sample10_50gen_anc.txt")) #invasion high (fixed for one ancestry!)
  plot_anc(make_dip("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/inv_high_neg1/sample10_50gen_anc.txt")) #invasion high = sliver (selected for light green, but mostly overwhelmed by dark green demographics)
  plot_anc(make_dip("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/inv_high_pos1/sample10_50gen_anc.txt")) #invasion high - selection for invaders (Dark green)
  plot_anc(make_dip("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/inv_high_pos1/sample10_250gen_anc.txt")) #invasion high (fixed for one ancestry!) #this looks light but I think it should be dark green (!) that they all fix for
  plot_anc(make_dip("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/iso4_1/sample10_250gen_anc.txt"))  
  
  #2chrom simulations
  plot_anc(make_dip("2chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/iso4_1/sample10_250gen_anc_chr1.txt"))
  plot_anc(make_dip("2chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/iso4_1/sample10_250gen_anc_chr2.txt"))  
  plot_anc(make_dip("2chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/pairs_low_1/sample10_250gen_anc_chr1.txt"))
  plot_anc(make_dip("2chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/pairs_low_1/sample10_250gen_anc_chr2.txt"))  
  plot_anc(make_dip("2chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/panmixia_1/sample10_250gen_anc_chr1.txt"))
  plot_anc(make_dip("2chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/panmixia_1/sample10_250gen_anc_chr2.txt"))    
  plot_corr_matrix("2chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/iso4_1/sample10_250gen_anc_chr1.txt")
  plot_corr_matrix("2chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/iso4_1/sample10_250gen_anc_chr2.txt")  
  plot_corr_matrix("2chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/pairs_low_1/sample10_250gen_anc_chr1.txt")
  plot_corr_matrix("2chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/pairs_low_1/sample10_250gen_anc_chr2.txt")  
  plot_corr_matrix("2chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/panmixia_1/sample10_250gen_anc_chr1.txt")
  plot_corr_matrix("2chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/panmixia_1/sample10_250gen_anc_chr2.txt")    
  
  plot_corr_matrix("2chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/iso4_1/sample10_50gen_anc_chr2.txt")  
    
  plot_corr_matrix(input_file)
  plot_corr_matrix("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/step4_low_neutral/sample10_250gen_anc.txt", ind = T) #no clear pattern of shared drift in deviation from individual alpha
  plot_corr_matrix("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/step4_low_neutral/sample10_250gen_anc.txt", ind = F) #much clearer pattern of shared ancestry dev. from sample-wide alpha
  plot_corr_matrix("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/step4_high_neutral/sample10_250gen_anc.txt", ind = F) #ack!what's happening with some of the high mig. simulations
  plot_corr_matrix("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/iso4_neutral/sample10_250gen_anc.txt", ind = F) #nice, clear
  plot_corr_matrix("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/cont_mig_high_neutral/sample10_250gen_anc.txt", ind = F) #not working
  plot_corr_matrix("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/cont_mig_low_neutral/sample10_250gen_anc.txt") #pretty random - as expected
  plot_corr_matrix("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/panmixia_neutral/sample10_250gen_anc.txt", ind = T) #nice
  plot_corr_matrix("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/panmixia_neutral/sample10_250gen_anc.txt", ind = T) #nice
  plot_corr_matrix("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/inv_high_neutral/sample10_250gen_anc.txt", ind = F)#really not working
  plot_corr_matrix("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/inv_high_neutral/sample10_50gen_anc.txt")
  #works but not what I was expecting, can see higher covariance in pop 3 and between 3 and 4, but invasion starts in pop 3 with gene flow towards 6 so I expected more covariance between pops 5 & 6 than 3 & 4
  
  plot_corr_matrix("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/inv_high_neg1/sample10_50gen_anc.txt") #some not working
  plot_corr_matrix("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/pairs_high_neutral/sample10_250gen_anc.txt")#nice, though faint
  plot_corr_matrix("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/pairs_low_neutral/sample10_250gen_anc.txt") #more clear
  plot_corr_matrix("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/pairs_high_1/sample10_250gen_anc.txt") #overall pattern not affected much by selection
  plot_corr_matrix("1chrom/anc_haplotypes/perfect_ancestry_20000hap_10000mark/ms_sim_output/pairs_low_1/sample10_250gen_anc.txt") #pattern genomewide more affeted by selection in this sim
}

##################################################3
# cholesky decomposition in R
library(MASS)
test_Chol = function(){
  sigma_K = matrix(c(1,.2,.3, .2, 1, -.5, .3, -.5, 1), nrow = 3, ncol = 3)
  mu = c(1, 5, 2)
  mvn = mvrnorm(n=50000, # create some MVN data
                mu = mu, 
                Sigma = sigma_K)
  hist(mvn) #not normal - a mixture of correlated normals
  cor(mvn[,1], mvn[,2]) # empirical is close to expected .2
  chol_mvn = t(chol(sigma_K)) # take the transpose to get a lower triangular matrix
  chol_mvn %*% t(chol_mvn) == sigma_K #Sigma = LL'
  #inv_chol_mvn = solve(chol_mvn)
  inv_chol_mvn = calcInvL(sigma_K)
  z_mvn = apply(mvn, 1, function(i) inv_chol_mvn %*% (i - mu))
  cor(z_mvn[1,], z_mvn[2,]) # no correlation anymore, great
  cor(z_mvn[3,], z_mvn[2,])
  cor(z_mvn[3,], z_mvn[1,])
  hist(z_mvn) # ~N(0,1)
  hist(rnorm(50000*3, mean=0, sd=1), col="darkblue", add=TRUE)
  hist(z_mvn[1,])# ~ N(0,1) 
  # so it works well using the true sigma_K.
  # what if I estimate K from the simulated data?
  mu_est = colMeans(mvn) # mu_est = alpha instead of known mu
  sigma_K_est = calcK(ancFreqMatrix = t(mvn), 
                      alpha = mu_est)
  inv_chol_mvn_est = calcInvL(sigma_K_est)
  z_mvn_est = apply(mvn, 1, function(i) inv_chol_mvn_est %*% (i - mu_est))
  hist(z_mvn_est)
  # the estimates of mu_est and sigma_K_est and nearly spot on to the true simulated model values (many observations)
  
  # but what about a case where I'm constrained to stay within 0-1 ancestry frequencies
  # and what if I only get to observe ancestry for a small set of haplotypes (vs. the true population ancestry frequency as above)
  N = 10 # ten haplotypes per population
  #bin <- 
  }

ZAnc = function(ancFreq, invL, alpha){
  # ancestry is a vector of raw ancestry frequency for each population
  # invL is inverse of (the lower triangular) Cholesky decomposition of the ancestry variance-covariance matrix
  # alpha is the population mean ancestries (a vector), based off of some neutral region or genomewide avg.
  
  # throws an error if the dimensions are not correct (all should = total # of populations)
  if (!(length(alpha) == dim(invL)[1] && 
        length(alpha) == length(ancFreq))){
    stop("Oops! dimensions for ancFreq, invL and alpha must match")
  }
  sum(zAnc(ancFreq, invL, alpha))
}

# same transformation as ZAnc but doesn't sum across populations
zAnc = function(ancFreq, invL, alpha){
  # throws an error if the dimensions are not correct (all should = total # of populations)
  if (!(length(alpha) == dim(invL)[1] && 
        length(alpha) == length(ancFreq))){
    stop("Oops! dimensions for ancFreq, invL and alpha must match")
  }
  invL %*% (ancFreq - alpha)
} #inv_chol_mvn %*% (i - mu))

zBeta = function(ancFreq, envWeights, invL, alpha){
  # calculates the association between transformed ancestry frequencies
  # and transformed environmental weights (e.g. c(1,0,1,0) if discrete) or c(altitude1, altitude2) if continuous
  # based on the inverse cholesky of the population frequencies
  # and where the environment is a selective environment (or apriori assumed relative strength of selection within the set of environments)
  if (!(length(alpha) == dim(invL)[1] && 
        length(alpha) == length(ancFreq) &&
        length(alpha) == length(envWeights))){
    stop("Oops! dimensions for envWeights, ancFreq, invL and alpha must match")
  }
  zAnc = invL %*% (ancFreq - alpha)
  zEnv = invL %*% (envWeights - mean(envWeights))
  # how well is the transformed ancestry predicted by its corresponding transformed environment?
  fit = lm(zAnc ~ zEnv)
  # summary(fit)
  fit$coefficients["zEnv"] # output slope of simple regression
}

# zBeta2 returns more information about the linear model proposed: 
# (the intercept, the estimated slope, sum of squared residuals, adjusted R-squared, and p-value)
# it also does NOT mean center the environment
zBeta2 = function(ancFreq, envWeights, invL, alpha){
  # calculates the association between transformed ancestry frequencies
  # and transformed environmental weights (e.g. c(1,0,1,0) if discrete) or c(altitude1, altitude2) if continuous
  # based on the inverse cholesky of the population frequencies
  # and where the environment is a selective environment (or apriori assumed relative strength of selection within the set of environments)
  if (!(length(alpha) == dim(invL)[1] && 
        length(alpha) == length(ancFreq) &&
        length(alpha) == length(envWeights))){
    stop("Oops! dimensions for envWeights, ancFreq, invL and alpha must match")
  }
  zAnc = invL %*% (ancFreq - alpha)
  zEnv = invL %*% (envWeights - mean(envWeights)) # centers environmental weights
  # how well is the transformed ancestry predicted by its corresponding transformed environment?
  fit = lm(zAnc ~ zEnv)
  fit.summary = summary(fit)
  # summary(fit)
  return(c(intercept = fit$coefficients[1],
           slope = fit$coefficients[2],
           sum_sq_res = sum(fit$residuals^2),
           r_squared = fit.summary$r.squared, # coefficient of determination
           pval_intercept = fit.summary$coefficients[1,4],
           pval_slope = ifelse(is.na(fit$coefficients[2]), NA, fit.summary$coefficients[2,4])))
}


# zBeta3 does NOT mean center the environment vector
# and it takes away the untransformed intercept -- this is the model I should use going forwards
# with option zInt =T it adds a transformed intercept = invL %*% vector of 1's
# note: transformed intercept should be used for environmental correlations, but not used for 1/0 vectors
# for selection present or not
zBeta3 = function(ancFreq, envWeights, invL, alpha, zInt = F){
  # calculates the association between transformed ancestry frequencies
  # and transformed environmental weights (e.g. c(1,0,1,0) if discrete) or c(altitude1, altitude2) if continuous
  # based on the inverse cholesky of the population frequencies
  # and where the environment is a selective environment (or apriori assumed relative strength of selection within the set of environments)
  if (!(length(alpha) == dim(invL)[1] && 
        length(alpha) == length(ancFreq) &&
        length(alpha) == length(envWeights))){
    stop("Oops! dimensions for envWeights, ancFreq, invL and alpha must match")
  }
  
  zAnc = invL %*% (ancFreq - alpha)
  zEnv = invL %*% envWeights
  
  # how well is the transformed ancestry predicted by its corresponding transformed environment?
  if (zInt){ # use a transformed intercept (i.e. ancestry freqs can have an environmental association and be higher or lower than the mean)
    zInt = invL %*% rep(1, length(envWeights)) # make intercept
    fit = lm(zAnc ~ 0 + zEnv + zInt) # transformed intercept but no untransformed intercept
    fit.summary = summary(fit)
    output = c(fit$coefficients[1:2],
               sum_sq_res = sum(fit$residuals^2),
               r_squared = fit.summary$r.squared, # coefficient of determination
               pval_zEnv = fit.summary$coefficients[1, 4],
               pval_zInt = tryCatch(fit.summary$coefficients[2, 4], error = function(e){
                 print("Warning: slope and intercept are confounded. try using no zInt intercept") 
                 return(NA)})) # if zEnv is something like c(1,1,1,1) R will only estimate one slope and throw an error
  } else{
    fit = lm(zAnc ~ 0 + zEnv) # no intercept
    fit.summary = summary(fit)
    output = c(fit$coefficients[1],
             zInt = NA,
             sum_sq_res = sum(fit$residuals^2),
             r_squared = fit.summary$r.squared, # coefficient of determination
             pval_zEnv = fit.summary$coefficients[1,4],
             pval_zInt = NA)
  }
  
  return(output)
}




ancBeta = function(ancFreq, envWeights, alpha){
  # calculates the association between untransformed ancestry frequencies
  # and untransformed environmental weights (e.g. c(1,0,1,0) if discrete) or c(altitude1, altitude2) if continuous
  # where the environment is a selective environment (or apriori assumed relative strength of selection within the set of environments)
  if (!(length(alpha) == length(ancFreq) &&
        length(alpha) == length(envWeights))){
    stop("Oops! dimensions for envWeights, ancFreq and alpha must match")
  }
  devAnc = (ancFreq - alpha)
  devEnv = (envWeights - mean(envWeights))
  # how well is the transformed ancestry predicted by its corresponding transformed environment?
  fit = lm(devAnc ~ devEnv)
  # summary(fit)
  fit$coefficients["devEnv"] # output slope of simple regression
}



calcProbZAnc = function(ZAnc, nPop, logProb = F){
  pnorm(ZAnc, mean = 0, sd = sqrt(nPop), 
        # variance = # of pops
        lower.tail = F, log.p = logProb)
}

# Poisson Binomial probability calculation (of = or more extreme - larger - value)
calcProbPoiBin = function(ancCount, nHap, alpha){ 
  # AncCount is a vector of counts w/ length = # pops
  # nHap is the number of haplotypes per pop, alpha is a vector of mean ancestry per pop
  # lower inclusive probability <= kk total counts
  lower = ppoibin(kk = sum(ancCount), pp = alpha, wts = rep(nHap, length(alpha)))
  # wts gives number of binomial draws (haplotypes) per alpha (population included)
  exact = dpoibin(kk = sum(ancCount), pp = alpha, wts = rep(nHap, length(alpha)))
  upper_equal = (1 - lower) + exact
  upper_equal
}

NoRot = function(K, ancFreq, alpha){# like ZAnc but with no rotation, 
  #just a non-rotating standardization (result = standard non-indep. normal ancestry freqs)
  sum((ancFreq - alpha)/diag(sqrt(K)))
}


simple_env_regression = function(ancFreq, envWeights){
  # simple regression between an environment (centered) and set of ancestry frequencies.
  if (!(length(envWeights) == length(ancFreq))){
    stop("Oops! dimensions for envWeights, ancFreq, invL and alpha must match")
  }
  fit = lm(ancFreq ~ envWeights)
  fit.summary = summary(fit)
  output = c(fit$coefficients[1:2],
    sum_sq_res = sum(fit$residuals^2),
    r_squared = fit.summary$r.squared, # coefficient of determination
    pval_Env = fit.summary$coefficients[1, 4],
    pval_Int = tryCatch(fit.summary$coefficients[2, 4], error = function(e){
      print("Warning: slope and intercept are confounded. try using no intercept") 
      return(NA)}))
  return(output)
}

  
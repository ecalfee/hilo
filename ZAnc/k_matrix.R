suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(poibin))

# functions
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

calcInvL = function(K){
  solve(t(chol(K)))  
  #calculate the inverse Cholesky. t(chol(K)) is a lower left triangular matrix such that LL' = K
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

# cholesky transformation of population ancestry frequencies
zAnc = function(ancFreq, invL, alpha){
  # throws an error if the dimensions are not correct (all should = total # of populations)
  if (!(length(alpha) == dim(invL)[1] && 
        length(alpha) == length(ancFreq))){
    stop("Oops! dimensions for ancFreq, invL and alpha must match")
  }
  invL %*% (ancFreq - alpha)
} #inv_chol_mvn %*% (i - mu))


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

# truncate to [0,1] range (e.g. from simulated mvn data)
truncate01 <- function(x){
  t <- x
  t[t < 0] <- 0
  t[t > 1] <- 1
  return(t)
}
# function to return log likelihood for MVN given data y
ll_mvn <- function(y, mu, detK, invK){ # here y is a vector of pop ancestry frequencies
  J = length(y) # number of data points
  ll = -1/2*log(detK) - J/2*log(2*pi) - 1/2*t(y - mu) %*% invK %*% (y - mu)
  return(ll)
}
# find maximum likelihood b under a model of environmental selection
# precalculate invK = solve(K)
ML_b <- function(y, alpha, invK, X){
  b = solve(t(X) %*% invK %*% X) %*% t(X) %*% invK %*% (y - alpha)
  return(b)
}

# AIC with small sample size correction
AICc <- function(ll, k, n){# ll is the log likelihood, k number of estimated parameters, and n the number of data points (here pops)
  2*k - 2*ll + 2*k*(k + 1)/(n - k - 1)
}


##################################################
# cholesky decomposition in R -- simple example

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
}

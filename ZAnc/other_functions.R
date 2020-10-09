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
  return(b[,]) # return as a scalar
}

# AIC with small sample size correction
AICc <- function(ll, k, n){# ll is the log likelihood, k number of estimated parameters, and n the number of data points (here pops)
  2*k - 2*ll + 2*k*(k + 1)/(n - k - 1)
}

# residual sum of squared errors
# if mu = alpha, this is the Anc_LK statistic
RSS <- function(ancFreq, mu, invK){
  rss = t(ancFreq - mu) %*% invK %*% (ancFreq - mu)
  return(rss[,]) # return as value, not 1x1 matrix
}

# fit the neutral model
fit_null <- function(anc, alpha, invK, detK){
  #invK = solve(K)
  #detK = det(K) # determinant of K matrix
  
  # estimate log likelihood under this MVN model
  logliks <- sapply(1:nrow(anc), function(i)
    ll_mvn(t(anc[i, ]), # make into a column vector
           mu = alpha,
           detK = detK,
           invK = invK))
  
  rss <- sapply(1:nrow(anc), function(i)
    RSS(ancFreq = t(anc[i, ]),
        mu = alpha,
        invK = invK))
  
  fits <- data.frame(
    RSS = rss,
    ll = logliks,
    k = 0,
    n = length(alpha)) %>%
    mutate(AICc = AICc(ll = ll, k = k, n = n))
  
  return(fits)
}

# fit a selection model with environment X:
fit_sel <- function(anc, alpha, detK, invK, X, b_names){
  # X is a matrix but may only have 1 column (e.g. X = [1,1,1..1]^T)
  # b_names can be a single value or vector, e.g. "b" or c("b0", "b1")
  if (!is.matrix(X)){
    stop("environment X must be a matrix")
  }
  #invK = solve(K)
  #detK = det(K)
  
  # estimate beta for each observed vector of population ancestries
  betas = apply(anc, 1, function(y)
    ML_b(y = y, alpha = alpha, invK = invK, X = X)) %>%
    as.matrix(.)
  if (dim(betas)[2] != length(b_names)) betas = t(betas) # fix because for 1D it outputs one way and 2D it outputs a diff way from apply
  #betas = apply(maize_anc, 1, function(y)
  #  ML_b(y = y, alpha = alpha, invK = invK, X = X)) %>%
  #  as.matrix(.)
  logliks <- sapply(1:nrow(anc), function(i)
    ll_mvn(t(anc[i, ]), # make into a column vector
           mu = alpha + X %*% betas[i, ], # use ML beta estimate calculate expected value
           detK = detK,
           invK = invK))
  
  rss <- sapply(1:nrow(anc), function(i)
    RSS(ancFreq = t(anc[i, ]),
        mu = alpha + X %*% betas[i, ],
        invK = invK))
  
  fits <- betas %>% # put everything together
    as.data.frame(.) %>%
    data.table::setnames(b_names) %>%
    mutate(RSS = rss,
           ll = logliks,
           k = length(b_names), # no error is estimated, free parameters = number of betas
           n = length(alpha),
           AICc = AICc(ll = ll, k = k, n = n))
  return(fits)
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

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

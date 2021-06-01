# functions to calculate false discovery rates
fdr_high <- function(a, data, sims){# high FDR for value a (expected/obs that are at least a)
  obs <- sum(data >= a)
  null <- sum(sims >= a)/length(sims)*length(data)
  null/obs
}
fdr_low <- function(a, data, sims){# low FDR for value a
  obs <- sum(data <= a)
  null <- sum(sims <= a)/length(sims)*length(data)
  null/obs
}

# Calculate 1% 5% and 10% FDRs
FDR_range = c(0.1, 0.05, 0.01)

# function to calculate high tail FDR
calc_FDR_high <- function(d, s, FDR_values, test_values){
  # calc (high) FDR for every value in test_values
  test_fdr_high <- sapply(test_values, function(x) 
    fdr_high(x, data = d, sims = s))
  # find which values match the FDR thresholds in FDR_values
  FDRs_high <- data.frame(threshold = sapply(FDR_values, 
                                             function(p) 
                                               min(test_values[test_fdr_high <= p], 
                                                   na.rm = T)),
                          FDR = FDR_values,
                          tail = "high",
                          stringsAsFactors = F)
  return(FDRs_high)
}

# function to calculate low tail FDR
calc_FDR_low <- function(d, s, FDR_values, test_values){
  # calc (low) FDR for every value in test_values
  test_fdr_low <- sapply(test_values, function(x) 
    fdr_low(x, data = d, sims = s))
  # find which values match the FDR thresholds in FDR_values
  FDRs_low <- data.frame(threshold = sapply(FDR_values, 
                                            function(p) 
                                              max(test_values[test_fdr_low <= p], 
                                                  na.rm = T)),
                         FDR = FDR_values,
                         tail = "low",
                         stringsAsFactors = F)
  return(FDRs_low)
}

# putting it together - get both high and low FDRs
calc_FDR <- function(d, s, FDR_values = FDR_range, 
                     test_values = seq(0, 1, by = .0001)){
  FDRs_low <- calc_FDR_high(d = d, s = s, 
                           FDR_values = FDR_values,
                           test_values = test_values)
  FDRs_high <- calc_FDR_low(d = d, s = s, 
                            FDR_values = FDR_values,
                            test_values = test_values)
  FDRs <- rbind(FDRs_high, FDRs_low)
  return(FDRs)
}

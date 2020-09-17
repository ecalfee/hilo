# function to fit simple linear model ancestry ~ environment
simple_env_regression = function(ancFreq, envWeights){
  if (!(length(envWeights) == length(ancFreq))){
    stop("Oops! dimensions for envWeights and ancFreq must match")
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
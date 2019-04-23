# script to plot ancestry frequency trajectories as described by 
# equations in Juric's supplement:
# analytical solution for admixture proportion at 
# time t after a single pulse
single_pulse <- function(t, a0, r, s){
  G <- (1-r)*(1-s)
  at <- a0*((G^t)*(1 - r - G) + r)/(1 - G)
  return(at)
}
# approximation for long t
single_pulse_inf_approx <- function(a0, r, s){
  a0*r/(1 - (1 - r)*(1 - s))
}
# approximation for long t and small r/s
single_pulse_inf_small_rs_approx <- function(a0, r, s){
  a0*r/(r + s)
}
# iterative solution one generation at a time
# for single pulse admixture event
single_pulse_iter <- function(t, a0, r, s){
  ts <- 0:t
  # store haplotype frequencies
  x <- integer(t+1) # neutral A ancestry, del minor A ancestry at linked selected site
  y <- integer(t+1) # neutral A ancestry, non-del major B ancestry at linked selected site
  # set initial conditions
  x[1] <- a0 # AA migrants from admixture pulse
  y[1] <- 0 # no AB ancestry haplotypes in generation 0
  for (i in 2:(t+1)){
    x[i] <- x[i-1]*(1-r)*(1-s)
    y[i] <- x[i-1]*r + y[i-1]
  }
  a = x + y # frequency of A ancestry at neutral locus
  return(data.frame(a=a, x=x, y=y, stringsAsFactors = F))
}

plot_single_pulse <- function(ts = 1:100, a0, r, s){
  plot(ts, single_pulse(ts, a0, r, s),
       type = "l", ylab = "a", xlab = "t",
       ylim = c(0, .5),
       main = "expected admixed ancestry proportion over time")
  abline(h = single_pulse_inf_approx(a0, r, s), col = "red")
  abline(h = single_pulse_inf_small_rs_approx(a0, r, s), col = "blue")
}
# plot a few examples:
# analytical solution vs. approximations:
plot_single_pulse(ts = 1:100, a0 = .2, r = .01, s=.01)
# iterative model matches analytical solution

plot(1:100, 
     single_pulse_iter(t = 100, a0 = .2, r = .01, s=.01)$a, 
     col = "green", type = "l",
     ylim = c(0, 0.5),
     main = "iterative admixture proportion over time")
lines(ts, 
      single_pulse_iter(t = 100, a0 = .2, r = .01, s=.01)$x, 
      col = "blue")
lines(ts, 
      single_pulse_iter(t = 100, a0 = .2, r = .01, s=.01)$y,
      col = "yellow")
lines(1:100, single_pulse(1:100, a0=.2, r=.01, s=.01), lty=2)
abline(h = single_pulse_inf_approx(a0=.2, r=.01, s=.01), col = "red", lty = 2)
abline(h = single_pulse_inf_small_rs_approx(a0=.2, r=.01, s=.01), col = "orange", lty = 2)



# continuous migration model
cont_mig <- function(t, a0 = .1, m, r, s, n = 0){
  G <- (1-r)*(1-s)
  x0 <- a0 # all migrants have AA haplotypes
  y0 <- 0 # no starting AB haplotypes
  xt <- x0*(G^t + m*(G^t - G)/(G - 1))
  yt <- x0*r*(1 + (G^t - G)/(G - 1) + 
                m*(G^t - t*G + t - 1)/(G - 1)^2)
  at <- xt + yt # total # A neutral alleles = AA + AB ancestry haplotypes
  return(at)
}

# iterative solution one generation at a time
# for continuous migration admixture
cont_mig_iter_orig <- function(t, a0, m, r, s, n = 0){
  ts <- 0:t
  # store haplotype frequencies
  x <- integer(t+1) # neutral A ancestry, del minor A ancestry at linked selected site
  y <- integer(t+1) # neutral A ancestry, non-del major B ancestry at linked selected site
  # set initial conditions
  x0 <- a0 # AA migrants from initial admixture pulse
  y0 <- 0 # no AB ancestry haplotypes in generation 0
  x[1] <- x0 
  y[1] <- y0
  for (i in 2:(t+1)){
    x[i] <- (x[i-1] + m*x0)*(1-r)*(1-s)
    y[i] <- (x[i-1] + m*x0)*r + y[i-1]
  }
  a = x + y # frequency of A ancestry at neutral locus
  return(data.frame(a=a, x=x, y=y, stringsAsFactors = F))
}

# plot
plot(0:10, cont_mig_iter_orig(t = 10, a0 = .1, m = .01, 
                               r = .1, s =.1, 
                               n = 0)$a,
     col = "green", type = "l", 
     main = "continuous migration model",
     xlab = "t", ylab = "a")
lines(0:10, cont_mig(t = 0:10, a0 = .1,
                     m = .01, 
                        r = .1, s = .1, 
                        n = 0), lty = 2)


# plot migration = a0 each generation (explodes)
plot(0:100, cont_mig_iter_orig(t = 100, a0 = .001, m = 1, 
                              r = .1, s =.1, 
                              n = 0)$a,
     col = "green", type = "l", 
     main = "continuous migration model",
     xlab = "t", ylab = "a")
lines(0:100, cont_mig(t = 0:100, a0 = .001,
                     m = 1, 
                     r = .1, s = .1, 
                     n = 0), lty = 2)


# different continuous migration model:
# migration is constant m each generation from minor ancestry
# and constant proportion n from major ancestry
# (note: m is not a proportion of a0 now, but replacement fraction)

cont_mig_iter_new <- function(t, m, r, s, n = 0){
  ts <- 0:t
  # store haplotype frequencies
  x <- integer(t+1) # neutral A ancestry, del minor A ancestry at linked selected site
  y <- integer(t+1) # neutral A ancestry, non-del major B ancestry at linked selected site
  # set initial conditions
  x[1] <- 0 # start with no A ancestry in AA
  y[1] <- 0 # or AB haplotypes
  for (i in 2:(t+1)){
    # 1 - m - n fraction don't get replaced from prior generation
    # and fraction m come in this generation with x haplotype (=AA)
    x[i] <- (x[i-1]*(1 - m - n) + m)*(1-r)*(1-s)
    y[i] <- (x[i-1]*(1 - m - n) + m)*r + y[i-1]*(1 - m - n)
  }
  a = x + y # frequency of A ancestry at neutral locus
  return(data.frame(a=a, x=x, y=y, stringsAsFactors = F))
}

plot(0:1000, cont_mig_iter_new(t = 1000, m = .01,
                              r = .1, s =.1,
                              n = .01)$a,
          col = "green", type = "l", 
          main = "new continuous migration model",
          xlab = "t", ylab = "a")
#abline()






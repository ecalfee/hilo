library(rethinking)
library(dplyr)
library(ggplot2)

# first make a smaller set of maize data
small_maize <- bind_cols(sites, maize_anc) %>%
  sample_n(., size = 1000, replace = F) %>%
  tidyr::gather(., "population", "mexicana_anc", maize_pops) %>%
  mutate(pop = as.integer(factor(population)))

m <- map( # quadratic approximation of the posterior MAP
  alist(
    mexicana_anc ~ dbeta2(prob = p, theta = theta), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- alpha[pop],
    alpha[pop] ~ dnorm(0, 10),
    #mu ~ dnorm(0, 5),
    theta ~ dunif(0, 10)
  ),
  data = small_maize#,
  #start = list(mu = 1, theta = 2, alpha = 0)
  )
precis(m, depth = 2)
pairs(m)
m.post <- sim(m, unique(small_maize[ , c("pop", "population")]), n = 1000)
m.post.alpha <- apply(m.post, 2, mean)
m.post.var <- apply(m.post, 2, var)
# how does the posterior estimated mean mexicana ancestry (alpha)
# compare to the data's mean, and the sample used to fit the model?
sample.alpha = group_by(small_maize, population) %>%
  summarise(alpha = mean(mexicana_anc)) %>%
  .$alpha
sample.var = group_by(small_maize, population) %>%
  summarise(var = var(mexicana_anc)) %>%
  .$var
plot(zAnc_maize$alpha, m.post.alpha, 
     main = "posterior model fit to alpha - pop mean")
points(zAnc_maize$alpha, sample.alpha, col = "purple")
legend("topleft", c("model", "sample"), 
       col = c("black", "purple"), pch = 1)
abline(0, 1, col = "blue")
# fits well for higher alpha and poorly for low alpha
# possibly because I need to give each alpha it's own variance.
# let's see:

# each pop has it's own variance
m1 <- map( # quadratic approximation of the posterior MAP
  alist(
    mexicana_anc ~ dbeta2(prob = p, theta = theta[pop]), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- alpha[pop],
    alpha[pop] ~ dnorm(0, 10),
    theta[pop] ~ dunif(0, 10)
  ),
  data = small_maize#,
)
precis(m1, depth = 2)
coef(m1)
pairs(m1)
m1.post <- sim(m1, unique(small_maize[ , c("pop", "population")]), n = 1000)
m1.post.alpha <- apply(m1.post, 2, mean)
m1.post.var <- apply(m1.post, 2, var)

# how does the posterior estimated mean mexicana ancestry (alpha)
# compare to the data's mean, and the sample used to fit the model?
plot(zAnc_maize$alpha, m1.post.alpha,
     main = "posterior model fit to alpha - pop var & mean")
points(zAnc_maize$alpha, sample.alpha, col = "purple")
legend("topleft", c("model", "sample"), 
       col = c("black", "purple"), pch = 1)
abline(0, 1, col = "blue")
# how well do variances fit?
plot(diag(zAnc_maize$K), m1.post.var,
     main = "posterior model fit to var - pop var & mean")
points(diag(zAnc_maize$K), sample.var, col = "purple")
legend("topleft", c("model", "sample"), 
       col = c("black", "purple"), pch = 1)
abline(0, 1, col = "blue")

# individual population distributions look really good
hist(small_maize$mexicana_anc[small_maize$population=="pop360"], col = "purple")
hist(m1.post[ , 1], add = T)
hist(rbeta2(1000, logistic(coef(m1)["alpha[1]"]), theta = coef(m1)["theta[1]"]),
     add = T)
hist(small_maize$mexicana_anc[small_maize$population=="pop361"], col = "purple")
hist(m1.post[ , 2], add = T)
hist(rbeta2(1000, logistic(coef(m1)["alpha[2]"]), theta = coef(m1)["theta[2]"]),
     add = T)

# but this doesn't capture correlations between populations, so I expect
# the distribution of the mean across populations to be off:
sample.global.means = small_maize %>%
  group_by(chr, pos) %>%
  summarise(mean_mexicana = mean(mexicana_anc)) %>%
  .$mean_mexicana
# global mean is very off
hist(sample.global.means, col = "purple")
hist(apply(m1.post, 1, mean), add = T)
# because it doesn't capture covariances:
plot(small_maize$mexicana_anc[small_maize$pop == 5],
     small_maize$mexicana_anc[small_maize$pop == 4])
abline(lm(small_maize$mexicana_anc[small_maize$pop == 5] ~
       small_maize$mexicana_anc[small_maize$pop == 4]), col = "blue")
cor(small_maize$mexicana_anc[small_maize$pop == 5],
     small_maize$mexicana_anc[small_maize$pop == 4])

# add covariances between admixed populations:
#sigmas <- c(sigma_a,sigma_b) # standard deviations
#Rho <- matrix( c(1,rho,rho,1) , nrow=2 ) # correlation matrix
# now matrix multiply to get covariance matrix
#Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)
#R ∼ LKJcorr(14) # correlation matrix

# start with just 2 populations: 4 & 5
small_maize_45 <- filter(small_maize, pop %in% c(4,5)) %>%
  dplyr::select(-pop) %>%
  tidyr::spread(., population, mexicana_anc)
m2 <- map( # quadratic approximation of the posterior MAP
  alist(
    mexicana_anc ~ dbeta2(prob = p, theta = theta[pop]), # dbeta2 is a reformulation of the 
    # beta distribution provided by the 'rethinking' package
    # where mean probability 'prob' = a/(a+b) and parameter 'theta' = a+b
    # are used instead of shape1 = a and shape2 = b
    logit(p) <- value[pop],
    value[pop] ~ dmvnorm(alpha[pop], , )
    alpha_pop365 ~ dnorm(0, 10),
    theta_pop365 ~ dunif(0, 10),
    Rho_pop ~ dlkjcorr(2),
    sigma_dept ~ dcauchy(0,2),
  ),
  data = small_maize_45
)
precis(m2, depth = 2)
coef(m2)
pairs(m2)
m2.post <- sim(m2, unique(small_maize_45[ , c("pop", "population")]), n = 1000)
m2.post.alpha <- apply(m2.post, 2, mean)
m2.post.var <- apply(m2.post, 2, var)
plot(m2.post[, 1], m2.post[, 2])
abline(lm(m2.post[, 2] ~
            m2.post[, 1]), col = "black")
points(small_maize$mexicana_anc[small_maize$pop == 5],
     small_maize$mexicana_anc[small_maize$pop == 4],
     col = "purple")
abline(lm(small_maize$mexicana_anc[small_maize$pop == 5] ~
            small_maize$mexicana_anc[small_maize$pop == 4]), col = "purple")



### MVN version
m3 <- map( # quadratic approximation of the posterior MAP
  alist(
    #pop363 ~ dnorm(alpha_pop363, theta_pop363),
    #c(pop363, pop365) ~ dmvnorm2(Mu = c(alpha_pop363, alpha_pop365), 
    #                        sigma = c(theta_pop363, theta_pop365),
    #                        Rho = diag(2)),
                            #Rho_pop),
   # c(a,b) ~ dmvnorm(c(alpha_pop363, alpha_pop365),
  #                              c(theta_pop363, theta_pop365) * diag(2) * c(theta_pop363, theta_pop365)),
    c(a, b) ~ dmvnorm2(Mu = c(alpha_pop363, alpha_pop365), 
                            sigma = c(theta_pop363, theta_pop365),
                            Rho = diag(2)),
    alpha_pop363 ~ dnorm(0, 10),
    alpha_pop365 ~ dnorm(0, 10),
    theta_pop363 ~ dunif(0, 10),
    theta_pop365 ~ dunif(0, 10)#,
    #Rho_pop ~ dlkjcorr(2) # weakly regularizing prior on correlation matrix
  ),
  data = small_maize_45
    )

m3 <- map( # quadratic approximation of the posterior MAP
  alist(
    mexicana_anc ~ dbeta2(prob = p, theta[pop]),
    logit(p) <- mu[pop],
    mu[pop] ~ dmvnorm2(c(alpha4, alpha5), 1, Rho),
    #alpha[pop] ~ dnorm(0, 10),
    alpha4 ~ dnorm(0, 10),
    alpha5 ~ dnorm(0, 10),
    theta[pop] ~ dunif(0, 10),
    Rho ~ dlkjcorr(2) # weakly regularizing prior on correlation matrix
  ),
  data = small_maize_45
)
precis(m3, depth = 2)
coef(m3)
pairs(m3)
m3.post <- sim(m3, unique(small_maize_45[ , c("pop", "population")]), n = 1000)
m3.post.alpha <- apply(m3.post, 2, mean)
m3.post.var <- apply(m3.post, 2, var)
plot(m3.post[, 1], m3.post[, 2])
abline(lm(m3.post[, 2] ~
            m3.post[, 1]), col = "black")
points(small_maize$mexicana_anc[small_maize$pop == 5],
       small_maize$mexicana_anc[small_maize$pop == 4],
       col = "purple")
abline(lm(small_maize$mexicana_anc[small_maize$pop == 5] ~
            small_maize$mexicana_anc[small_maize$pop == 4]), col = "purple")

# example
data("Kline2") # load data for total tools per island
data(islandsDistMatrix) # load data for distance between islands
Kline2 # small dataset of total tools per population
Kline2$society <- 1:10
mKline1 <- map( # quadratic approximation (instead of bayesian est.)
  alist(
    total_tools ~ dpois(lambda),
    log(lambda) <- a + g[society] + bp*logpop,
    g[society] ~ GPL2( Dmat , etasq , rhosq , 0.01 ),
    a ~ dnorm(0,10),
    bp ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1)
  ),
  data=list(
    total_tools=Kline2$total_tools,
    logpop=Kline2$logpop,
    society=Kline2$society,
    Dmat=islandsDistMatrix)
  )
mKline2 <- map2stan( # bayesian example implementation
  alist(
    total_tools ~ dpois(lambda),
    log(lambda) <- a + g[society] + bp*logpop,
    g[society] ~ GPL2( Dmat , etasq , rhosq , 0.01 ),
    a ~ dnorm(0,10),
    bp ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1)
  ),
  data=list(
    total_tools=Kline2$total_tools,
    logpop=Kline2$logpop,
    society=Kline2$society,
    Dmat=islandsDistMatrix),
  warmup=2000 , iter=1e4 , chains=4 )
mKline2.post <- sim(mKline2, n = 10)
# add a made-up second observation for each island
Kline3 <- Kline2 %>%
  mutate(total_tools = mKline2.post[ , 1]) %>%
  mutate(observation = 2) %>%
  bind_rows(., mutate(Kline2, observation = 1))
mKline3 <- map2stan( # bayesian example implementation
  alist(
    total_tools ~ dpois(lambda),
    log(lambda) <- a + g[society] + bp*logpop,
    g[society] ~ GPL2( Dmat , etasq , rhosq , 0.01 ),
    a ~ dnorm(0,10),
    bp ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1)
  ),
  data=list(
    total_tools=Kline3$total_tools,
    logpop=Kline3$logpop,
    society=Kline3$society,
    observation=Kline3$observation,
    Dmat=islandsDistMatrix),
  warmup=2000 , iter=1e4 , chains=4 )
precis(mKline3)
# Kline 4 creates a matrix of observations, where columns are societies and rows are obs.
Kline4 <- Kline2 %>%
  mutate(total_tools = mKline2.post[ , 1]) %>%
  mutate(observation = 2) %>%
  bind_rows(., mutate(Kline2, observation = 1)) %>%
  dplyr::select(culture, total_tools, observation) %>%
  tidyr::spread(., culture, total_tools) %>%
  dplyr::select(c("observation", colnames(islandsDistMatrix)))
mKline4 <- map2stan( # does not work
  alist(
    total_tools ~ dpois(lambda[observation]),
    log(lambda[observation]) <- a + g + bp*logpop,
    g ~ GPL2( Dmat , etasq , rhosq , 0.01 ),
    a ~ dnorm(0,10),
    bp ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1)
  ),
  data=list(
    total_tools=Kline4[ , 2:11],
    observation=Kline4$observation,
    Dmat=islandsDistMatrix),
  warmup=2000 , iter=1e4 , chains=4 )


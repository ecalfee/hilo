# this script runs a set of test cases using the empirical variance/covariance matrix for maize ancestry
# and example scenarios of selection (on/off) or correlation with environment (e.g. elevation)
# to load data needed to run (e.g. K matrix), first run the first part of zanc_selection_models.R

# test set: there are 14 populations of maize. I will create some test sets with 14 pops
# that have perfect environmental correlations
# each row is a made-up ancestry frequency observation; each column is a population 
test_case_names <- c("all_1", "last_pop_only", "sixth_pop_only", "first_pop_only", "linear_env", "first_one_out", "last_one_out", "linear_env+100", "linear_env*2", "linear_env*2+100")

# note that the ancestries I make up as test cases ignore any covariances between populations -- this isn't a perfect test,
# but still the correct model should be closer than any of the other proposed models.
test_anc <- rbind(rep(1, 14), # all selected for mexicana ancestry
                  c(zAnc_maize$alpha[1:13], 1), # only last population selected
                  c(zAnc_maize$alpha[1:5], 1, zAnc_maize$alpha[7:14]), # 6th population selected
                  c(1, zAnc_maize$alpha[2:14]), # last population selected
                  seq(0, 1, length.out = 14),
                  c(0, rep(1, 13)),
                  c(rep(1,13), 0),
                  seq(0, 1, length.out = 14) + 100,
                  seq(0, 1, length.out = 14)*2,
                  seq(0, 1, length.out = 14)*2 + 100) # linear association with simple linear environmental gradient: e.g. seq(100, 200, length.out = 14)
# could also test false positive cases: e.g. one outlier but no real other association with

# I'll use my test_anc cases with the empirical K matrix:
test_env <- rbind(rep(1, 14), # all selected for mexicana ancestry
                  c(zAnc_maize$alpha[1:13], 1), # only last population selected
                  c(zAnc_maize$alpha[1:5], 1, zAnc_maize$alpha[7:14]), # 6th population selected
                  c(1, zAnc_maize$alpha[2:14]), # last population selected
                  seq(0, 1, length.out = 14)) # linear environment

# transformed test environments:
apply(test_env, 1, function(e) zAnc_maize$InvL %*% e)
apply(test_env, 1, function(e) zAnc_maize$InvL %*% (e - mean(e)))

zTest2 <- lapply(1:nrow(test_env), function(e) # no mean centering of environment
  apply(test_anc, 1, function(l) # calculate slope for each locus
    zBeta2(ancFreq = l,
           envWeights = test_env[e, ], 
           invL = zAnc_maize$InvL, 
           alpha = zAnc_maize$alpha)))

zTest3_noInt <- lapply(1:nrow(test_env), function(e) # takes away untransformed intercept too
  apply(test_anc, 1, function(l) # calculate slope for each locus
    zBeta3(ancFreq = l,
           envWeights = test_env[e, ], 
           invL = zAnc_maize$InvL, 
           alpha = zAnc_maize$alpha,
           zInt = F)))
a2 <- do.call(rbind,
              lapply(1:length(zTest2), function(l)
                data.frame(t(zTest2[[l]])) %>%
                  mutate(test_model = test_case_names) %>%
                  mutate(true_model = test_case_names[l])))
a2 %>%
  ggplot(aes(x = test_model, y = r_squared, color = true_model, shape = (true_model == test_model))) +
  geom_point() + 
  ggtitle("mean centered environment test cases") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
a2 %>%
  gather(., "model_output", "value", c("intercept..Intercept.", "slope.zEnv", "r_squared")) %>%
  ggplot(aes(x = test_model, y = value, color = true_model, shape = (true_model == test_model))) +
  geom_point() + 
  ggtitle("mean centered environment test cases") +
  facet_wrap(~model_output) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/test_cases_environment_ZAnc_model_values_centered_env.png",
       height = 10, width = 10, device = "png", units = "in")

a2 %>%
  gather(., "model_output", "value", c("sum_sq_res", "pval_intercept", "pval_slope.zEnv")) %>%
  mutate(., log10_value = log10(value)) %>%
  ggplot(aes(x = test_model, y = log10_value, color = true_model, shape = (true_model == test_model))) +
  geom_point() + 
  ggtitle("mean centered environment test cases") +
  facet_wrap(~model_output) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/test_cases_environment_ZAnc_fit_values_centered_env.png",
       height = 10, width = 10, device = "png", units = "in")


a3_noInt <- do.call(rbind,
              lapply(1:length(zTest3_noInt), function(l)
                data.frame(t(zTest3_noInt[[l]])) %>%
                  mutate(test_model = test_case_names) %>%
                  mutate(true_model = test_case_names[l])))
a3_noInt %>%
  ggplot(aes(x = test_model, y = r_squared, color = true_model, shape = (true_model == test_model))) +
  geom_point() + 
  ggtitle("not mean centered environment test cases") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
a3_noInt %>%
  gather(., "model_output", "value", c("zEnv", "r_squared")) %>%
  mutate(., log10_value = log10(value)) %>%
  ggplot(aes(x = true_model, y = log10_value, color = test_model, shape = (true_model == test_model))) +
  geom_point() + 
  ggtitle("not mean centered environment test cases") +
  facet_wrap(~model_output) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/test_cases_environment_ZAnc_model_values.png",
       height = 10, width = 10, device = "png", units = "in")

a3_noInt %>%
  gather(., "model_output", "value", c("sum_sq_res", "pval_zEnv")) %>%
  mutate(., log10_value = log10(value)) %>%
  ggplot(aes(x = true_model, y = log10_value, color = test_model, shape = (true_model == test_model))) +
  geom_point() + 
  ggtitle("not mean centered environment test cases") +
  facet_wrap(~model_output) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/test_cases_environment_ZAnc_fit_values.png",
       height = 10, width = 10, device = "png", units = "in")

# add in an intercept of 1's vectors for the environmental correlation case
# b/c I don't want the environmental effect to be constrained to the non-centered environment
# but I won't test the c(1,1,1,1) model with this
zTest3_int <- lapply(1:nrow(test_env), function(e)
  apply(test_anc[ , ], 1, function(l) # calculate slope for each locus
    zBeta3(ancFreq = l,
               envWeights = test_env[e, ], 
               invL = zAnc_maize$InvL, 
               alpha = zAnc_maize$alpha,
               zInt = T)))


# hmm..problem is that r-squared doesn't have the same meaning if there is an intercept or not
# (all models are being compared to free intercept model, not zTz now)
a3_int <- do.call(rbind,
                  lapply(1:length(zTest3_int), function(l)
                    data.frame(t(zTest3_int[[l]])) %>%
                      mutate(test_model = test_case_names) %>%
                      mutate(true_model = test_case_names[l])))
a3_int %>%
  gather(., "model_output", "value", c("zEnv", "zInt", "r_squared")) %>%
  #mutate(., log10_value = log10(value)) %>%
  ggplot(aes(x = true_model, y = value, color = test_model, shape = (true_model == test_model))) +
  geom_point() + 
  ggtitle("not mean centered environment test cases + transformed intercept") +
  facet_wrap(~model_output, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/test_cases_environment_ZAnc_model_values_3int.png",
       height = 10, width = 10, device = "png", units = "in")

a3_int %>% # here last pop only and first pop only are basically equivalent models because of the rotated intercept adding a second degree of freedom
  gather(., "model_output", "value", c("sum_sq_res", "pval_zEnv")) %>%
  mutate(., log10_value = log10(value)) %>%
  ggplot(aes(x = true_model, y = log10_value, color = test_model, shape = (true_model == test_model))) +
  geom_point() + 
  ggtitle("not mean centered environment test cases + intercept") +
  facet_wrap(~model_output) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/test_cases_environment_ZAnc_fit_values_3int.png",
       height = 10, width = 10, device = "png", units = "in")

# now plot each model with it's appropriate intercept (zElev) or no intercept (101001 selection/no-sel models):
# the true model always has the lowest sum_sq_error and linear permutations of the true model are equivalent (for environment)
zTest3 <- lapply(1:nrow(test_env), function(e)
  apply(test_anc[ , ], 1, function(l) # calculate slope for each locus
    zBeta3(ancFreq = l,
           envWeights = test_env[e, ], 
           invL = zAnc_maize$InvL, 
           alpha = zAnc_maize$alpha,
           zInt = c(F, F, F, F, T)[e]))) #only use intercept for linear environment model

a3 <- do.call(rbind,
                  lapply(1:length(zTest3), function(l)
                    data.frame(t(zTest3[[l]])) %>%
                      mutate(test_model = test_case_names) %>%
                      mutate(true_model = test_case_names[l])))
a3 %>%
  gather(., "model_output", "value", c("zEnv", "zInt", "r_squared")) %>%
  #filter(., value <= 5) %>%
  ggplot(aes(x = true_model, y = value, color = test_model, shape = (true_model == test_model))) +
  geom_point() + 
  ggtitle("correct model: zIntercept only for linear environment.") +
  facet_wrap(~model_output, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/test_cases_environment_ZAnc_model_values_correct3.png",
       height = 10, width = 10, device = "png", units = "in")
a3 %>% # good that sum of sq residuals pulls out the correct model (or linear equivalent in the case of linear environment model) 
  gather(., "model_output", "value", c("sum_sq_res", "pval_zEnv", "pval_zInt")) %>%
  mutate(., log10_value = log10(value)) %>%
  ggplot(aes(x = true_model, y = log10_value, color = test_model, shape = (true_model == test_model))) +
  geom_point() + 
  ggtitle("correct model: zIntercept only for linear environment") +
  facet_wrap(~model_output) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/test_cases_environment_ZAnc_fit_values_correct3.png",
       height = 10, width = 10, device = "png", units = "in")


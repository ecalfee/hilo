# in this script I plot sensitivity to the prior
# by comparing mean ancestry from the local ancestry inference
# to (1) the prior and (2) global ancestry estimates from NGSAdmix
MAIN_DIR = "../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm"

# TO DO: Instead, I should be plotting an individual's prior and ancestry output, with size = coverage,
# to get a sense of what coverage is necessary to reduce bias. Also really need to put in .99 and .01 prior comparisons
# GRAHAM's suggestion was also to mask sites per individual where they do not have any high-confidence call (maybe low information there)
# so, for example, I could look at just the sites in the genome where one of their posterior's is > .85


# meta data on populations
pops <- unique(read.table("../data/pass1_ids.txt", stringsAsFactors = F, 
                          header = T, sep = "\t")[ , c("popN", "zea", "LOCALITY")])
  
#tot_coverage_incl <- 

# 'true' global ancestry estimated from NGSAdmix
global <- read.table(paste0(MAIN_DIR, "/input/globalAdmixtureByPopN.txt"),
                     stringsAsFactors = F, header = F)
colnames(global) <- c("popN", "alpha_maize", "alpha_mex")

# equal ancestry ratios as prior
prior_equal <- do.call(rbind, lapply(global$popN, function(i)
  read.table(paste0(MAIN_DIR, "/output_equal_prior_noBoot/anc/pop", i, ".alpha.ind"),
             stringsAsFactors = F, header = T) %>%
    dplyr::mutate(popN = i))) %>%
  dplyr::group_by(popN) %>%
  dplyr::summarise(alpha_ancestryHMM = mean(alpha), n = n()) %>%
  dplyr::mutate(prior = "equal_ancestry")

# very extreme high mexicana ancestry as prior (.999 vs .001)
prior_very_high_mex <- do.call(rbind, lapply(global$popN, function(i)
  read.table(paste0(MAIN_DIR, "/output_very_high_mex_prior_noBoot/anc/pop", i, ".alpha.ind"),
             stringsAsFactors = F, header = T) %>%
    dplyr::mutate(popN = i))) %>%
  dplyr::group_by(popN) %>%
  dplyr::summarise(alpha_ancestryHMM = mean(alpha), n = n()) %>%
  dplyr::mutate(prior = "very_high_mex")

# very extreme high maize ancestry as prior (.999 vs .001)
prior_very_high_maize <- do.call(rbind, lapply(global$popN, function(i)
  read.table(paste0(MAIN_DIR, "/output_very_high_maize_prior_noBoot/anc/pop", i, ".alpha.ind"),
             stringsAsFactors = F, header = T) %>%
    dplyr::mutate(popN = i))) %>%
  dplyr::group_by(popN) %>%
  dplyr::summarise(alpha_ancestryHMM = mean(alpha), n = n()) %>%
  dplyr::mutate(prior = "very_high_maize")

# high maize ancestry as prior (.99 vs .01)
prior_high_maize <- do.call(rbind, lapply(global$popN, function(i)
  read.table(paste0(MAIN_DIR, "/output_high_maize_prior_noBoot/anc/pop", i, ".alpha.ind"),
             stringsAsFactors = F, header = T) %>%
    dplyr::mutate(popN = i))) %>%
  dplyr::group_by(popN) %>%
  dplyr::summarise(alpha_ancestryHMM = mean(alpha), n = n()) %>%
  dplyr::mutate(prior = "high_maize")

# high mex ancestry as prior (.99 vs .01)
prior_high_mex <- do.call(rbind, lapply(global$popN, function(i)
  read.table(paste0(MAIN_DIR, "/output_high_mex_prior_noBoot/anc/pop", i, ".alpha.ind"),
             stringsAsFactors = F, header = T) %>%
    dplyr::mutate(popN = i))) %>%
  dplyr::group_by(popN) %>%
  dplyr::summarise(alpha_ancestryHMM = mean(alpha), n = n()) %>%
  dplyr::mutate(prior = "high_mex")


#sapply(global$popN, function(i) mean(read.table(paste0(MAIN_DIR, "/output_equal_prior_noBoot/anc/pop", i, ".alpha.ind"),
#                             stringsAsFactors = F, header = T)$alpha))
# using 50/50 prior
prior_NGSadmix <- do.call(rbind, lapply(global[global$alpha_maize > 0 & global$alpha_mex > 0, "popN"], function(i)
  read.table(paste0(MAIN_DIR, "/output_noBoot/anc/pop", i, ".alpha.ind"),
             stringsAsFactors = F, header = T) %>%
    dplyr::mutate(popN = i))) %>%
  dplyr::group_by(popN) %>%
  dplyr::summarise(alpha_ancestryHMM = mean(alpha), n = n()) %>%
  dplyr::mutate(prior = "NGSAdmix_estimates")

d <- bind_rows(left_join(prior_NGSadmix, global, by = "popN"), left_join(prior_equal, global, by = "popN"),
               left_join(prior_very_high_maize, global, by = "popN"), left_join(prior_very_high_mex, global, by = "popN"),
               left_join(prior_high_maize, global, by = "popN"), left_join(prior_high_mex, global, by = "popN")) %>%
  left_join(., pops, by = "popN")

p1 <- d %>%
  ggplot(., aes(x = alpha_mex, y = alpha_ancestryHMM, color = prior)) +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_point() +
  ggtitle("Effect of prior on mean ancestry estimates from HMM")
plot(p1)
ggsave("plots/alpha_est_NGSadmix_vs_ancestry_hmm_compare_priors.png", 
       plot = p1,
       height = 10, width = 10, units = "in", device = "png")


# plot alpha NGSAdmix vs. alpha from ancestry_hmm just for NGSAdmix prior
p2 <- d %>%
  filter(., prior == "NGSAdmix_estimates") %>%
  ggplot(., aes(x = alpha_mex, y = alpha_ancestryHMM)) +
  geom_point(aes(col = LOCALITY)) +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  ggtitle("genome-wide mexicana ancestry: NGSadmix vs. ancestry_hmm w/ NGSadmix prior")
plot(p2)
ggsave("plots/alpha_est_NGSadmix_vs_ancestry_hmm_w_NGSadmix_prior.png", 
       plot = p2,
       height = 5, width = 5, units = "in", device = "png")

# plot alpha NGSAdmix vs. alpha from ancestry_hmm just for equal prior
p3 <- d %>%
  filter(., prior == "equal_ancestry") %>%
  ggplot(., aes(x = alpha_mex, y = alpha_ancestryHMM)) +
  geom_point(aes(col = LOCALITY)) +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  ggtitle("genome-wide mexicana ancestry: NGSadmix vs. ancestry_hmm w/ 50/50 prior")
plot(p3)
ggsave("plots/alpha_est_NGSadmix_vs_ancestry_hmm_w_equal_ancestry_prior.png", 
       plot = p3,
       height = 5, width = 5, units = "in", device = "png")
  
  
  
  
# in this script I analyze ancestry covariances between individuals
# and populations to make inferrences about demography
# and selection from genome-wide patterns.

# later I need to exclude large inversions (!)

# get ancestry, K matrices and useful functions
source("ZAnc_statistic.R")
source("melt_k_matrix.R")

# ancestry covariance plots
# population data
# plot population K as a correlation matrix
png(paste("../plots/K_corrplot_combined.png"), # saves plot as pin in ../plots/
    height = 7, width = 7, units = "in", res = 150)
corrplot(cov2cor(with_inv$pop$all$K), method="shade", main = "mex/maize pop anc corr.")
dev.off()

# plot population K raw covariance matrix
png(paste("../plots/K_covplot_combined.png"), # saves plot as pin in ../plots/
    height = 7, width = 7, units = "in", res = 150)
#heatmap(all$K, Colv = NA, Rowv = NA, main = "mex/maize pop anc cov.",
#        col= colorRampPalette(brewer.pal(8, "Blues"))(25))
col_pop_location=location_colors[my_group=as.numeric(as.factor(pops$LOCALITY))]
heatmap(with_inv$pop$all$K, Colv = NA, Rowv = NA, main = "mex/maize pop anc cov.",
        col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
        RowSideColors = col_pop_location,
        scale = "none")
dev.off()

# individual data
# plot individual K, where individuals are grouped (for visualization) by pop,
# as a correlation matrix
png(paste("../plots/K_corrplot_ind_combined.png"), # saves plot as pin in ../plots/
    height = 21, width = 21, units = "in", res = 150)
corrplot(cov2cor(with_inv$ind$all$K), method="shade", main = "mex/maize ind anc corr.")
dev.off()

# plot individual K as a raw covariance matrix
png(paste("../plots/K_covplot_ind_combined.png"), # saves plot as pin in ../plots/
    height = 21, width = 21, units = "in", res = 150)
col_ind_location=location_colors[my_group = 
                                   as.numeric(as.factor(included_inds$LOCALITY))]
heatmap(with_inv$ind$all$K, 
        Colv = NA, Rowv = NA, 
        main = "mex/maize ind anc cov.",
        col = colorRampPalette(brewer.pal(8, "Blues"))(25),
        RowSideColors = col_ind_location,
        scale = "none")
dev.off()

# I don't have a good idea yet how much of the ancestry covariance across all maize has to do with selection
# at the same loci for mexicana introgression.
# comparing K matrices for low and high recombining regions of the genome
# will help answer this question of selection vs. demography

# A first step is to get recombination rates for every thinned position
# and plot mean ancestry by recombination rate 
# and visualize K by high/low recomb. where alpha is a genomewide mean 
# and visualize K by high/low recomb. where alpha is a recombination bin mean

# low vs. high recomb. population level K 
png(paste("../plots/K_covplot_pop_combined_low_recomb_rate.png"), # saves plot as pin in ../plots/
    height = 7, width = 7, units = "in", res = 150)
heatmap(with_inv$pop$low_r$K, Colv = NA, Rowv = NA, main = "low recomb - pop anc cov.",
        col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
        RowSideColors = col_pop_location,
        scale = "none")
dev.off()

png(paste("../plots/K_covplot_pop_combined_high_recomb_rate.png"), # saves plot as pin in ../plots/
    height = 7, width = 7, units = "in", res = 150)
heatmap(with_inv$pop$high_r$K, Colv = NA, Rowv = NA, main = "high recomb - pop anc cov.",
        col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
        RowSideColors = col_pop_location,
        scale = "none")
dev.off()


# low vs. high recomb. ind. level K
png(paste("../plots/K_covplot_ind_combined_low_recomb_rate.png"), # saves plot as pin in ../plots/
    height = 14, width = 14, units = "in", res = 150)
heatmap(with_inv$ind$low_r$K, Colv = NA, Rowv = NA, main = "low recomb - ind anc cov.",
        col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
        RowSideColors = col_ind_location,
        scale = "none")
dev.off()

png(paste("../plots/K_covplot_ind_combined_high_recomb_rate.png"), # saves plot as pin in ../plots/
    height = 14, width = 14, units = "in", res = 150)
heatmap(with_inv$ind$high_r$K, Colv = NA, Rowv = NA, main = "high recomb - ind anc cov.",
        col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
        RowSideColors = col_ind_location,
        scale = "none")
dev.off()

# correlation matrices low vs. high recomb.
# population level
png(paste("../plots/K_corrplot_pop_combined_low_recomb_rate.png"),
    height = 7, width = 7, units = "in", res = 150)
corrplot(cov2cor(with_inv$pop$low_r$K), method="shade", main = "Low recomb. pop anc corr.")
dev.off()
png(paste("../plots/K_corrplot_pop_combined_high_recomb_rate.png"),
    height = 7, width = 7, units = "in", res = 150)
corrplot(cov2cor(with_inv$pop$high_r$K), method="shade", main = "High recomb. pop anc corr.")
dev.off()



# individual level
png(paste("../plots/K_corrplot_ind_combined_low_recomb_rate.png"),
    height = 14, width = 14, units = "in", res = 150)
corrplot(cov2cor(with_inv$ind$low_r$K), method="shade", main = "Low recomb. ind anc corr.")
dev.off()
png(paste("../plots/K_corrplot_ind_combined_high_recomb_rate.png"),
    height = 14, width = 14, units = "in", res = 150)
corrplot(cov2cor(with_inv$ind$high_r$K), method="shade", main = "High recomb. ind anc corr.")
dev.off()

# do pairs of maize-mex populations at the same location 
# covary more with each other than random pairs?
# which populations have mex-maize pairs?
pop_pairs = pops$LOCALITY[duplicated(pops$LOCALITY)]
pops$paired = pops$LOCALITY %in% pop_pairs
pop_pairs_SHORT = unique(pops$LOC_SHORT[pops$LOCALITY %in% pop_pairs])

# I can plot the covariances from the K matrices 
# as box plots to compare covariances
# for specific pairs of pops

# melt or flatten k matrices to plot
melted_ks <- list()
melted_ks$pop <- do.call(rbind,
                         # for covariance data with and without the inversion
                         # melt covariance matrices into flat dataframes
                         # and merge the results
                         lapply(c("with_inv", "no_inv"), 
                                function(k) {
                                  ks <- get(k)
                                  results <- do.call(rbind, 
                                                     lapply(names(ks$pop), # all high_r low_r
                                                            function(i) kmelt(ks$pop[[i]]$K, 
                                                                              haps_per_pop = pops$N_haps, 
                                                                              correlation = F, regions = i)))
                                  results$inv <- k
                                  return(results)
                                }
                         ))
melted_ks$ind <- do.call(rbind,
                         # for covariance data with and without the inversion
                         # melt covariance matrices into flat dataframes
                         # and merge the results
                         lapply(c("with_inv", "no_inv"), 
                                function(k) {
                                  ks <- get(k)
                                  results <- do.call(rbind, 
                                                     lapply(names(ks$ind), # all high_r low_r
                                                            function(i) kmelt_ind(ks$ind[[i]]$K, 
                                                                              correlation = F, regions = i)))
                                  results$inv <- k
                                  return(results)
                                }
                         ))
#melted_ks$pop <- do.call(rbind,
#                     lapply(names(with_inv$pop), # all high_r low_r
#                            function(i) kmelt(with_inv$pop[[i]]$K, 
#                                              haps_per_pop = pops$N_haps, 
#                                              correlation = F, regions = i)))
#melted_ks$ind <- do.call(rbind,
#                         lapply(names(with_inv$pop), # all high_r low_r
#                                function(i) kmelt_ind(with_inv$ind[[i]]$K,
#                                                  correlation = F, regions = i)))


# I also want to look at correlations not covariances
# because the variance differs between low and high r regions
# because mean ancestry is closer to 50/50 in high r regions.
melted_corrs <- list()
# same as above but always correlation = T for the kmelt function
melted_corrs$pop <- do.call(rbind,
                         lapply(c("with_inv", "no_inv"), 
                                function(k) {
                                  ks <- get(k)
                                  results <- do.call(rbind, 
                                                     lapply(names(ks$pop), # all high_r low_r
                                                            function(i) kmelt(ks$pop[[i]]$K, 
                                                                              haps_per_pop = pops$N_haps, 
                                                                              correlation = T, regions = i)))
                                  results$inv <- k
                                  return(results)
                                }
                         ))
melted_corrs$ind <- do.call(rbind,
                         lapply(c("with_inv", "no_inv"), 
                                function(k) {
                                  ks <- get(k)
                                  results <- do.call(rbind, 
                                                     lapply(names(ks$ind), # all high_r low_r
                                                            function(i) kmelt_ind(ks$ind[[i]]$K, 
                                                                                  correlation = T, regions = i)))
                                  results$inv <- k
                                  return(results)
                                }
                         ))


# are populations that simply have pairs different?
# doesn't seem to make a difference
melted_ks$pop %>%
  filter(., mex_vs_maize) %>%
  ggplot(., aes(y = covariance, x = have_pairs,
                fill = regions)) + 
  geom_boxplot() +
  facet_wrap(~inv) +
  ggtitle("covariance mex-maize for pops that do and don't have same location pairs")
ggsave(filename = "../plots/anc_cov_mex-maize_pops_no_pairs_location.png",
       device = "png", 
       height = 3, width = 8, 
       units = "in")

# show patterns for covarience between maize-mex pairs
melted_ks$pop %>%
  filter(., have_pairs & mex_vs_maize) %>%
  ggplot(., aes(y = covariance, x = same_location,
                fill = regions)) + 
  geom_boxplot() +
  facet_wrap(~inv) +
  ggtitle("covariance mex-maize pairs at same and diff locations")
ggsave(filename = "../plots/anc_cov_mex-maize_pops_same_diff_location.png",
       device = "png", 
       height = 3, width = 8, 
       units = "in")
for (i in c("all", "low_r", "high_r")){
  print(anova(lm(covariance ~ same_location,  # not sig. diff
                 data = filter(melted_ks$pop, inv == "no_inv" & have_pairs & mex_vs_maize & regions == i ))))
}
# but are high and low r different?
# only marginally sig. with inv and
# not at all when excluding inversions
lapply(c("no_inv", "with_inv"), function(i){
  print(i)
  anova(lm(covariance ~ same_location*regions,
                 data = filter(melted_ks$pop, inv == i, have_pairs & mex_vs_maize & regions != "all")))
})


# show patterns of covariance for maize-maize pairs
melted_ks$pop %>%
  filter(., (!mex_vs_maize) & zea_SHORT_1 == "maiz") %>%
  ggplot(., aes(y = covariance, x = same_location,
                fill = regions)) + 
  facet_wrap(~inv) +
  geom_boxplot() +
  ggtitle("covariance maize-maize pairs at same and diff locations")
ggsave(filename = "../plots/anc_cov_maize-maize_pops_same_diff_location.png",
       device = "png", 
       height = 3, width = 8, 
       units = "in")
anova(lm(covariance ~ same_location,  # more variance w/in pop than covariance between pops
         data = filter(melted_ks$pop, inv == "no_inv" & (!mex_vs_maize) & zea_SHORT_1 == "maiz" & regions == "all")))


# show patterns of covariance for mex-mex pairs
melted_ks$pop %>%
  filter(., (!mex_vs_maize) & zea_SHORT_1 == "mexi") %>%
  ggplot(., aes(y = covariance, x = same_location,
                fill = regions)) + 
  facet_wrap(~inv) +
  geom_boxplot() +
  ggtitle("covariance mex-mex pairs at same and diff locations")
ggsave(filename = "../plots/anc_cov_mex-mex_pops_same_diff_location.png",
       device = "png", 
       height = 3, width = 8, 
       units = "in")

# now I want to look at individual by individual covariance
# to confirm that within-population individuals have more 
# ancestry covariance than between populations

# plot covariances ind
melted_ks$ind %>%
  filter(., have_pairs & mex_vs_maize) %>%
  ggplot(., aes(y = covariance, x = same_location,
                fill = regions)) + 
  facet_wrap(~inv) +
  geom_boxplot() +
  ggtitle("covariance mex-maize ind's at same and diff locations")
ggsave(filename = "../plots/anc_cov_mex-maize_same_diff_location.png",
       device = "png", 
       height = 3, width = 6, 
       units = "in")

# maize- mex pairs - overall low covariance across zea groups
# same location has little effect but high vs. low r matters
for (i in c("with_inv", "no_inv")){
  print(i)
  print("covariance across locations")
  print(anova(lm(covariance ~ same_location,  # no diff
           data = filter(melted_ks$ind, inv == i & have_pairs & mex_vs_maize & regions == "all"))))
  print("covariance mexicana")
  print(anova(lm(covariance ~ regions*same_location,  # mex cov
           data = filter(melted_ks$ind, inv == i & mex_vs_maize & regions != "all"))))
  print("correlation mexicana")
  print(anova(lm(covariance ~ regions*same_location,  # mex corr
           data = filter(melted_corrs$ind, inv == i & mex_vs_maize & regions != "all"))))
  }


# Given that we don't expect inbreeding, 
# I don't think I'm properly correcting for sample size N
# for self comparisons (variance instead of covariance)
# within individual vs. between individuals same pop
melted_ks$ind %>%
  filter(., (!mex_vs_maize) & zea_1 == "maize" & same_location) %>%
  ggplot(., aes(y = covariance, x = same_ind,
                fill = regions)) + 
  geom_boxplot() +
  facet_wrap(~inv) +
  ggtitle("w/in pop cov. w/in & between maize ind's")

# same w/in vs. between individuals but with mexicana
melted_ks$ind %>%
  filter(., (!mex_vs_maize) & zea_1 == "mexicana" & same_location) %>%
  ggplot(., aes(y = covariance, x = same_ind,
                fill = regions)) + 
  geom_boxplot() +
  facet_wrap(~inv) +
  ggtitle("w/in pop cov. w/in & between mex ind's")


# now excluding within-individual comparisons:
# show patterns of covariance for maize-maize individuals
# from the same and different populations
melted_ks$ind %>%
  filter(., (!mex_vs_maize) & zea_1 == "maize" & !same_ind) %>%
  ggplot(., aes(y = covariance, x = same_location,
                fill = regions)) + 
  geom_boxplot() +
  facet_wrap(~inv) +
  ggtitle("covariance maize-maize ind's at same and diff locations")
ggsave(filename = "../plots/anc_cov_maize-maize_same_diff_location.png",
       device = "png", 
       height = 3, width = 8, 
       units = "in")

# more covariance w/in pops than between for all, low and high r regions
for (i in c("all_ind", "low_r_ind", "high_r_ind")){
  print(anova(lm(covariance ~ same_location,  
                 data = filter(melted_ks_ind, (!mex_vs_maize) & zea_1 == "maize" & !same_ind & regions == i))))
}
# also low r have higher covariance than high r regions:
# but OPPOSITE pattern if we exclude the inversion
# and no sig. interaction between recomb. and within vs. between pops
anova(lm(covariance ~ same_location*regions,  
         data = filter(melted_ks$ind, 
                       (!mex_vs_maize) & zea_1 == "maize" & 
                         !same_ind & regions != "all")))

# show patterns of covariance for mex-mex individuals
melted_ks$ind %>%
  filter(., (!mex_vs_maize) & zea_1 == "mexicana" & !same_ind) %>%
  ggplot(., aes(y = covariance, x = same_location,
                fill = regions)) + 
  geom_boxplot() +
  facet_wrap(~inv) +
  ggtitle("covariance mex-mex ind's at same and diff locations")
ggsave(filename = "../plots/anc_cov_mex-mex_same_diff_location.png",
       device = "png", 
       height = 3, width = 8, 
       units = "in")

for (i in c("all_ind", "low_r_ind", "high_r_ind")){
  print(anova(lm(covariance ~ same_location,  # more covariance w/in pops than between
                 data = filter(melted_ks_ind, (!mex_vs_maize) & zea_1 == "mexicana" & !same_ind & regions == i))))
}
# in mexicana (not maize) there is a sig. (but very small) interaction
# between recomb. rate and within-vs-between pops comparisons
# also effect of being in same pop is a lot stronger
anova(lm(covariance ~ same_location*regions,  
         data = filter(melted_ks$ind, 
                       (!mex_vs_maize) & zea_1 == "mexicana" & 
                         !same_ind & regions != "all")))

# correlations:
# mex-maize
melted_corrs$ind %>% # effect of low vs. high r only (more covariance, and especially corr, at low r -- is this due to the inversion?)
  filter(., have_pairs & mex_vs_maize) %>%
  ggplot(., aes(y = covariance, x = same_location,
                fill = regions)) + 
  geom_boxplot() +
  facet_wrap(~inv) +
  ggtitle("correlation mex-maize ind's at same and diff locations")
ggsave(filename = "../plots/anc_corr_mex-maize_same_diff_location.png",
       device = "png", 
       height = 3, width = 8, 
       units = "in")

# maize-maize
melted_corrs$ind %>% # effect of low vs. high r only (more covariance, and especially corr, at low r -- is this due to the inversion?)
  filter(., !mex_vs_maize & zea_1 == "maize" & !same_ind) %>%
  ggplot(., aes(y = covariance, x = same_location,
                fill = regions)) + 
  geom_boxplot() +
  facet_wrap(~inv) +
  ggtitle("correlation maize-maize ind's at same and diff locations")
ggsave(filename = "../plots/anc_corr_maize-maize_same_diff_location.png",
       device = "png", 
       height = 3, width = 8, 
       units = "in")


# mex-mex individuals
melted_corrs$ind %>% # effect of low vs. high r only (more covariance, and especially corr, at low r -- is this due to the inversion?)
  filter(., !mex_vs_maize & zea_1 == "mexicana" & !same_ind) %>%
  ggplot(., aes(y = covariance, x = same_location,
                fill = regions)) + 
  geom_boxplot() +
  facet_wrap(~inv) +
  ggtitle("correlation mex-mex ind's at same and diff locations")
ggsave(filename = "../plots/anc_corr_mex-mex_same_diff_location.png",
       device = "png", 
       height = 3, width = 8, 
       units = "in")

# maize correlations
for (i in c("all_ind", "low_r_ind", "high_r_ind")){
  print(anova(lm(covariance ~ same_location,  # more correlation between ind's w/in same pop than diff pop
                 data = filter(melted_corrs_ind, (!mex_vs_maize) & zea_1 == "maize" & !same_ind & regions == i))))
}

# mexicana correlations
for (i in c("all_ind", "low_r_ind", "high_r_ind")){
  print(anova(lm(covariance ~ same_location,  # more correlation between ind's w/in same pop than diff pop
                 data = filter(melted_corrs_ind, (!mex_vs_maize) & zea_1 == "mexicana" & !same_ind & regions == i))))
}

# low recombining regions have more covariance than high r regions
# and more correlations between diff maize ind's. No interaction with same of diff pop:
anova(lm(covariance ~ regions*same_location,  # maize cov
         data = filter(melted_ks$ind, (!mex_vs_maize) & zea_1 == "maize" & !same_ind & regions != "all")))
anova(lm(covariance ~ regions*same_location,  # maize corr
         data = filter(melted_corrs$ind, (!mex_vs_maize) & zea_1 == "maize" & !same_ind & regions != "all")))
# in mexicana there is an interaction with same location.
anova(lm(covariance ~ regions*same_location,  # mex cov
         data = filter(melted_ks$ind, (!mex_vs_maize) & zea_1 == "mexicana" & !same_ind & regions != "all")))
anova(lm(covariance ~ regions*same_location,  # mex corr
         data = filter(melted_corrs$ind, (!mex_vs_maize) & zea_1 == "mexicana" & !same_ind & regions != "all")))

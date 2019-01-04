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
corrplot(cov2cor(all$K), method="shade", main = "mex/maize pop anc corr.")
dev.off()

# plot population K raw covariance matrix
png(paste("../plots/K_covplot_combined.png"), # saves plot as pin in ../plots/
    height = 7, width = 7, units = "in", res = 150)
#heatmap(all$K, Colv = NA, Rowv = NA, main = "mex/maize pop anc cov.",
#        col= colorRampPalette(brewer.pal(8, "Blues"))(25))
col_pop_location=location_colors[my_group=as.numeric(as.factor(pops$LOCALITY))]
heatmap(all$K, Colv = NA, Rowv = NA, main = "mex/maize pop anc cov.",
        col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
        RowSideColors = col_pop_location,
        scale = "none")
dev.off()

# individual data
# plot individual K, where individuals are grouped (for visualization) by pop,
# as a correlation matrix
png(paste("../plots/K_corrplot_ind_combined.png"), # saves plot as pin in ../plots/
    height = 21, width = 21, units = "in", res = 150)
corrplot(cov2cor(all_ind$K), method="shade", main = "mex/maize ind anc corr.")
dev.off()

# plot individual K as a raw covariance matrix
png(paste("../plots/K_covplot_ind_combined.png"), # saves plot as pin in ../plots/
    height = 21, width = 21, units = "in", res = 150)
col_ind_location=location_colors[my_group = 
                                   as.numeric(as.factor(included_inds$LOCALITY))]
heatmap(all_ind$K, 
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
heatmap(low_r$K, Colv = NA, Rowv = NA, main = "low recomb - pop anc cov.",
        col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
        RowSideColors = col_pop_location,
        scale = "none")
dev.off()

png(paste("../plots/K_covplot_pop_combined_high_recomb_rate.png"), # saves plot as pin in ../plots/
    height = 7, width = 7, units = "in", res = 150)
heatmap(high_r$K, Colv = NA, Rowv = NA, main = "high recomb - pop anc cov.",
        col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
        RowSideColors = col_pop_location,
        scale = "none")
dev.off()


# low vs. high recomb. ind. level K
png(paste("../plots/K_covplot_ind_combined_low_recomb_rate.png"), # saves plot as pin in ../plots/
    height = 14, width = 14, units = "in", res = 150)
heatmap(low_r_ind$K, Colv = NA, Rowv = NA, main = "low recomb - ind anc cov.",
        col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
        RowSideColors = col_ind_location,
        scale = "none")
dev.off()

png(paste("../plots/K_covplot_ind_combined_high_recomb_rate.png"), # saves plot as pin in ../plots/
    height = 14, width = 14, units = "in", res = 150)
heatmap(high_r_ind$K, Colv = NA, Rowv = NA, main = "high recomb - ind anc cov.",
        col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
        RowSideColors = col_ind_location,
        scale = "none")
dev.off()

# correlation matrices low vs. high recomb.
# population level
png(paste("../plots/K_corrplot_pop_combined_low_recomb_rate.png"),
    height = 7, width = 7, units = "in", res = 150)
corrplot(cov2cor(low_r$K), method="shade", main = "Low recomb. pop anc corr.")
dev.off()
png(paste("../plots/K_corrplot_pop_combined_high_recomb_rate.png"),
    height = 7, width = 7, units = "in", res = 150)
corrplot(cov2cor(high_r$K), method="shade", main = "High recomb. pop anc corr.")
dev.off()



# individual level
png(paste("../plots/K_corrplot_ind_combined_low_recomb_rate.png"),
    height = 14, width = 14, units = "in", res = 150)
corrplot(cov2cor(low_r_ind$K), method="shade", main = "Low recomb. ind anc corr.")
dev.off()
png(paste("../plots/K_corrplot_ind_combined_high_recomb_rate.png"),
    height = 14, width = 14, units = "in", res = 150)
corrplot(cov2cor(high_r_ind$K), method="shade", main = "High recomb. ind anc corr.")
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
melted_ks <- do.call(rbind,
                     lapply(c("all", "low_r", "high_r"),
                            function(i) kmelt(get(i)$K, haps_per_pop = pops$N_haps, cov2corr = F, regions = i)))
# show patterns for covarience between maize-mex pairs
melted_ks %>%
  filter(., have_pairs & mex_vs_maize) %>%
  ggplot(., aes(y = covariance, x = same_location,
                fill = regions)) + 
  geom_boxplot() +
  ggtitle("covariance mex-maize pairs at same and diff locations")
for (i in c("all", "low_r", "high_r")){
  print(anova(lm(covariance ~ same_location,  # not sig. diff
                 data = filter(melted_ks, have_pairs & mex_vs_maize & regions == i))))
}


# show patterns of covariance for maize-maize pairs
melted_ks %>%
  filter(., (!mex_vs_maize) & zea_SHORT_1 == "maiz") %>%
  ggplot(., aes(y = covariance, x = same_location,
                fill = regions)) + 
  geom_boxplot() +
  ggtitle("covariance maize-maize pairs at same and diff locations")
anova(lm(covariance ~ same_location,  # more variance w/in pop than covariance between pops
         data = filter(melted_ks, (!mex_vs_maize) & zea_SHORT_1 == "maiz" & regions == "all")))


# show patterns of covariance for mex-mex pairs
melted_ks %>%
  filter(., (!mex_vs_maize) & zea_SHORT_1 == "mexi") %>%
  ggplot(., aes(y = covariance, x = same_location,
                fill = regions)) + 
  geom_boxplot() +
  ggtitle("covariance mex-mex pairs at same and diff locations")

# now I want to look at individual by individual covariance
# to confirm that within-population individuals have more 
# ancestry covariance than between populations

# covariances ind
melted_ks_ind <- do.call(rbind,
                         lapply(c("all_ind", "low_r_ind", "high_r_ind"),
                                function(i) kmelt_ind(get(i)$K, cov2cor = F, regions = i)))

melted_ks_ind %>%
  filter(., have_pairs & mex_vs_maize) %>%
  ggplot(., aes(y = covariance, x = same_location,
                fill = regions)) + 
  geom_boxplot() +
  ggtitle("covariance mex-maize ind's at same and diff locations")
anova(lm(covariance ~ same_location,  # no diff
         data = filter(melted_ks_ind, have_pairs & mex_vs_maize & regions == "all_ind")))
melted_corrs_ind %>% # effect of low vs. high r only (more covariance, and especially corr, at low r -- is this due to the inversion?)
  filter(., have_pairs & mex_vs_maize) %>%
  ggplot(., aes(y = covariance, x = same_location,
                fill = regions)) + 
  geom_boxplot() +
  ggtitle("correlation mex-maize ind's at same and diff locations")
# maize- mex pairs - overall low covariance across zea groups
# same location has little effect but high vs. low r matters
anova(lm(covariance ~ regions*same_location,  # mex corr
         data = filter(melted_ks_ind, mex_vs_maize & regions != "all")))
anova(lm(covariance ~ regions*same_location,  # mex corr
         data = filter(melted_corrs_ind, mex_vs_maize & regions != "all")))



# within individual vs. between individuals same pop
melted_ks_ind %>%
  filter(., (!mex_vs_maize) & zea_1 == "maize" & same_location) %>%
  ggplot(., aes(y = covariance, x = same_ind,
                fill = regions)) + 
  geom_boxplot() +
  ggtitle("w/in pop cov. w/in & between maize ind's")

# same w/in vs. between individuals but with mexicana
melted_ks_ind %>%
  filter(., (!mex_vs_maize) & zea_1 == "mexicana" & same_location) %>%
  ggplot(., aes(y = covariance, x = same_ind,
                fill = regions)) + 
  geom_boxplot() +
  ggtitle("w/in pop cov. w/in & between mex ind's")


# now excluding within-individual comparisons:
# show patterns of covariance for maize-maize individuals
# from the same and different populations
melted_ks_ind %>%
  filter(., (!mex_vs_maize) & zea_1 == "maize" & !same_ind) %>%
  ggplot(., aes(y = covariance, x = same_location,
                fill = regions)) + 
  geom_boxplot() +
  ggtitle("covariance maize-maize ind's at same and diff locations")
# more covariance w/in pops than between for all, low and high r regions
for (i in c("all_ind", "low_r_ind", "high_r_ind")){
  print(anova(lm(covariance ~ same_location,  
                 data = filter(melted_ks_ind, (!mex_vs_maize) & zea_1 == "maize" & !same_ind & regions == i))))
}
# also low r have higher covariance than high r regions:
# but no sig. interaction between recomb. and within vs. between pops
anova(lm(covariance ~ same_location*regions,  
         data = filter(melted_ks_ind, 
                       (!mex_vs_maize) & zea_1 == "maize" & 
                         !same_ind & regions != "all")))

# show patterns of covariance for mex-mex individuals
melted_ks_ind %>%
  filter(., (!mex_vs_maize) & zea_1 == "mexicana" & !same_ind) %>%
  ggplot(., aes(y = covariance, x = same_location,
                fill = regions)) + 
  geom_boxplot() +
  ggtitle("covariance mex-mex ind's at same and diff locations")
for (i in c("all_ind", "low_r_ind", "high_r_ind")){
  print(anova(lm(covariance ~ same_location,  # more covariance w/in pops than between
                 data = filter(melted_ks_ind, (!mex_vs_maize) & zea_1 == "mexicana" & !same_ind & regions == i))))
}
# in mexicana (not maize) there is a sig. (but very small) interaction
# between recomb. rate and within-vs-between pops comparisons
# also effect of being in same pop is a lot stronger
anova(lm(covariance ~ same_location*regions,  
         data = filter(melted_ks_ind, 
                       (!mex_vs_maize) & zea_1 == "mexicana" & 
                         !same_ind & regions != "all")))



# I also want to look at correlations not covariances
# because the variance differs between low and high r regions
# because mean ancestry is closer to 50/50 in high r regions.
# correlations:
melted_corrs_ind <- do.call(rbind,
                            lapply(c("all_ind", "low_r_ind", "high_r_ind"),
                                   function(i) kmelt_ind(get(i)$K, cov2cor = T, regions = i)))
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
         data = filter(melted_ks_ind, (!mex_vs_maize) & zea_1 == "maize" & !same_ind & regions != "all")))
anova(lm(covariance ~ regions*same_location,  # maize corr
         data = filter(melted_corrs_ind, (!mex_vs_maize) & zea_1 == "maize" & !same_ind & regions != "all")))
# in mexicana there is an interaction with same location.
anova(lm(covariance ~ regions*same_location,  # mex cov
         data = filter(melted_ks_ind, (!mex_vs_maize) & zea_1 == "mexicana" & !same_ind & regions != "all")))
anova(lm(covariance ~ regions*same_location,  # mex corr
         data = filter(melted_corrs_ind, (!mex_vs_maize) & zea_1 == "mexicana" & !same_ind & regions != "all")))


# pairs_k is a k matrix with covariances between mexicana pops (rows)
# and maize pops (cols), where diagonal is same location,
# and only pops with a pair are included in the matrix
for (i in c("all", "low_r", "high_r")){
  k = get(i)
  pairs_k_unordered = k$K[pops$zea == "mexicana" & pops$paired, 
                          pops$zea == "maize" & pops$paired]
  k$pairs_k = pairs_k_unordered[order(rownames(pairs_k_unordered)), 
                                order(colnames(pairs_k_unordered))]
  # make dataframe with columns 'paired', 'unpaired', 'nonpair'
  cov_pairs = diag(k$pairs_k)
  cov_nonpairs = k$pairs_k[col(k$pairs_k) != row(k$pairs_k)]
  boxplot(cov_nonpairs, cov_pairs,
          main = paste0(i, " mex-maize anc. cov."),
          xlab = "diff vs. same loc")
}



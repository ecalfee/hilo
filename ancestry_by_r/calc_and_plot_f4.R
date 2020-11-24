library(dplyr)
library(tidyr)
library(ggplot2)
library(boot)
# this script calculated f4-ratio to get an admixture proportion estimate
# for a sympatric population across recombination rate and coding density quintiles
# working directory "~/Documents/gitErin/hilo")
# note: W384 has no f4 data b/c it's a small ~10kb window
# with no tripsacum reads mapped. See last line of depthSample -- hilo$ tail -n 1 ancestry_by_r/results/f4/sympatric_maize/W384.depthSample

# load variables from Snakefile
sympatric_pop = snakemake@params[["sympatric_pop"]]
# sympatric_pop = "sympatric_maize"
# sympatric_pop = "sympatric_mexicana"
# zea is maize or mexicana. Only used for colors on plots. maize pops have index numbers > 100
zea = ifelse(sympatric_pop == "sympatric_maize" || as.integer(substr(sympatric_pop, 4, 7)) > 100, "maize", "mexicana") 

f4_num_file = snakemake@input[["f4_num"]]
# f4_num_file = paste0("ancestry_by_r/results/f4/", sympatric_pop, ".f4")
f4_denom_file = snakemake@input[["f4_denom"]]
# f4_denom_file = paste0("ancestry_by_r/results/f4/allopatric_maize.f4")
windows_file = snakemake@input[["windows"]]
# windows_file = "ancestry_by_r/results/map_pos_1cM_windows.txt"
colors_file = snakemake@input[["colors"]]
# colors_file = "colors.R"
png_r5 = snakemake@output[["png_r5"]]
# png_r5 = paste0("ancestry_by_r/plots/f4_allo_pop22_symp_", sympatric_pop, "_byr5.png")
png_cd5 = snakemake@output[["png_cd5"]]
# png_cd5 = paste0("ancestry_by_r/plots/f4_allo_pop22_symp_", sympatric_pop, "_bycd5.png")
png_r5_no_inv4m = snakemake@output[["png_r5_no_inv4m"]]
# png_r5_no_inv4m = paste0("ancestry_by_r/plots/f4_allo_pop22_symp_", sympatric_pop, "_byr5_noinv4m.png")
png_cd5_no_inv4m = snakemake@output[["png_cd5_no_inv4m"]]
# png_cd5_no_inv4m = paste0("ancestry_by_r/plots/f4_allo_pop22_symp_", sympatric_pop, "_bycd5_noinv4m.png")
png_f4_num_denom = snakemake@output[["png_f4_num_denom"]]
# png_f4_num_denom = paste0("ancestry_by_r/plots/f4_allo_pop22_symp_", sympatric_pop, "_num_denom.png")
n_boot = snakemake@params[["n_boot"]]
# n_boot = 1000
inv_file = snakemake@input[["inv_file"]]
# inv_file = "data/refMaize/inversions/knownInv_v4_coord.txt"

source(colors_file)

# population permutation used to calculate D-statistic:
# (ABBA-BABA)/(ABBA+BABA) in phylogeny (((parv, x), allo_mex), trip)
#d_symp <- read.table(paste0("ancestry_by_r/results/f4/", sympatric_pop, ".2.f4"),
#                     stringsAsFactors = F, sep = "\t", header = T)
# d-stat across all sites in genome:
#sum(d_symp$Numer)/sum(d_symp$Denom)


# population permutation used to calculate f4-ratio:
# f4(parv, allo_mex; X, trip) = E[(x_parv - x_allo_mex)*(x_X - x_trip)] in phylogeny (((parv, x), allo_mex), trip)
# where X is a sympatric population (numerator)
# or X is non-admixed allopatric maize (denominator)

# f4 ratio numerator
f4_num <- read.table(f4_num_file, stringsAsFactors = F, sep = "\t", header = T)
# f4 ratio denominator allopatric maize
f4_denom <- read.table(f4_denom_file, stringsAsFactors = F, sep = "\t", header = T)
# angsd outputs sum of (x_parv - x_allo_mex)*(x_X - x_trip)
# so I divide by number of sites to get the expectation
alpha = (sum(f4_num$Numer)/sum(f4_num$numSites))/(sum(f4_denom$Numer)/sum(f4_denom$numSites))
# genomewide the alpha estimate for % maize is ~ 58% in sympatric maize

# look at f4 across recombination rates:
# get windows
winds <- read.table(windows_file, header = T, stringsAsFactors = F, sep = "\t")
# get inversion coordinates
inv <- read.table(inv_file, header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("inv", "chr", "start", "end", "length"))

# which windows overlap inv4m?
winds$inv4m = winds$chr == inv$chr[inv$inv == "inv4m"] &
  winds$start < inv$end[inv$inv == "inv4m"] &
  winds$end > inv$start[inv$inv == "inv4m"]

convert_2_Mb_per_cM <- function(bin){
  start = substr(bin, 1, 1)
  end = substr(bin, nchar(bin), nchar(bin))
  middle = substr(bin, 2, nchar(bin) - 1) %>%
    strsplit(., split = ",") %>%
    unlist(.) %>%
    as.numeric(.)
  range = round(1/middle, 4) # 1cM*1/(cM/Mb) = Mb/cM, rounded to 4 decimal places
  return(paste0(start, range[2], ",", range[1], end))
}
# test:
# convert_2_Mb_per_cM("(2,5.8]")

# combine data
d_f4 <- group_by(f4_num, window) %>%
  summarise(num_stat = sum(Numer),
            num_sites = sum(numSites)) %>%
  full_join(., group_by(f4_denom, window) %>%
              summarise(denom_stat = sum(Numer),
                        denom_sites = sum(numSites)),
            by = "window") %>%
  left_join(winds, ., by = "window") %>%
  # report quintiles for # bp in a cM window (analogous to cd5, which is coding bp/cM)
  # instead or recombination rate, cM/Mb 
  dplyr::mutate(quintile_Mbp5 = 6 - quintile_r5,
                bin_Mbp5 = sapply(bin_r5, function(x) convert_2_Mb_per_cM(x))
  ) %>%
  ungroup()

# by recombination rate (inverse of bp density)
f4_by_r <- d_f4 %>%
  group_by(., bin_Mbp5, quintile_Mbp5) %>%
  filter(!is.na(num_sites)) %>% # eliminates 2 windows w/out data
  summarise(f4_num = sum(num_stat)/sum(num_sites),
            f4_denom = sum(denom_stat)/sum(denom_sites),
            f4_alpha = f4_num/f4_denom,
            f4_mex = 1 - f4_alpha, # alpha measured proportion maize, so 1-alpha is proportion mexicana
            n_sites_tot = (sum(num_sites) + sum(denom_sites))/2,
            n_windows = length(unique(window))) %>%
  arrange(quintile_Mbp5) %>%
  ungroup()
#View(f4_by_r)

# by coding density
f4_by_cd <- d_f4 %>%
  group_by(., bin_cd5, quintile_cd5) %>%
  filter(!is.na(num_sites)) %>% # eliminates 2 windows w/out data
  summarise(f4_num = sum(num_stat)/sum(num_sites),
            f4_denom = sum(denom_stat)/sum(denom_sites),
            f4_alpha = f4_num/f4_denom,
            f4_mex = 1 - f4_alpha, # alpha measured proportion maize, so 1-alpha is proportion mexicana
            n_sites_tot = (sum(num_sites) + sum(denom_sites))/2,
            n_windows = length(unique(window))) %>%
  arrange(quintile_cd5) %>%
  ungroup()
#View(f4_by_cd)
#View(f4_by_r)

# with(f4_by_r, cor(x = quintile_Mbp5, y = f4_mex, method = "spearman"))

# bootstrap
# function to calc stat
calc_cor_r <- function(d, i) { # d is data, i is indices
  b = d[i, ] %>%
    group_by(., bin_Mbp5, quintile_Mbp5) %>%
    summarise(f4_num = sum(num_stat)/sum(num_sites),
              f4_denom = sum(denom_stat)/sum(denom_sites),
              f4_alpha = f4_num/f4_denom,
              f4_mex = 1 - f4_alpha,
              n = n()) %>%
    arrange(quintile_Mbp5)
  return(c(spearman = with(b, cor(x = quintile_Mbp5, 
                                  y = f4_mex, 
                                  method = "spearman")),
           # also get alpha estimates for each quintile in the bootstrap sample
           unlist(pivot_wider(b, names_from = quintile_Mbp5, 
                              names_prefix = "Mbp", 
                              values_from = f4_mex, 
                              id_cols = quintile_Mbp5)),
           # and return number of windows per quintile in bootstrap sample
           unlist(pivot_wider(b, names_from = quintile_Mbp5, 
                              names_prefix = "n", 
                              values_from = n, 
                              id_cols = quintile_Mbp5))))
}

# test:
#filter(d_f4, !is.na(num_sites)) %>% calc_cor_r(., 1:nrow(.)) # all windows
#filter(d_f4, !is.na(num_sites) & !inv4m) %>%
#  calc_cor_r(., 1:nrow(.)) # all windows not overlapping inversion
#filter(d_f4, !is.na(num_sites)) %>%
#  calc_cor_r(., sample(1:nrow(.), replace = T))

# bootstrap across all windows
boot_r <- boot(data = filter(d_f4, !is.na(num_sites)), 
               statistic = calc_cor_r,
               sim = "ordinary",
               stype = "i",
               R = n_boot)

# bootstrap confidence interval for spearman's rank correlation:
boot.ci(boot_r, index = 1, type = c("norm", "basic", "perc"))

# bootstrap confidence interval for number of low recombination bins:
#boot.ci(boot_r, index = 7, type = c("norm", "basic", "perc"))

# bootstrap across all windows, but 
# exclude windows overlapping inv4m
boot_r_no_inv4m <- boot(data = filter(d_f4, !is.na(num_sites) & !inv4m), 
                        statistic = calc_cor_r,
                        sim = "ordinary",
                        stype = "i",
                        R = n_boot)

# bootstrap confidence interval for spearman's rank correlation:
#boot.ci(boot_r_no_inv4m, index = 1, type = c("norm", "basic", "perc"))

# bootstrap confidence interval for number of low recombination bins:
#boot.ci(boot_r_no_inv4m, index = 11, type = c("norm", "basic", "perc"))


# bootstrap sampling within each quintile separately:
# (produces a very similar answer)
# function to calc stat
calc_cor_r5 <- function(d, i) { # d is data, i is indices
  all = 1:nrow(d)
  j = ifelse(i == all, i, # original sample (not bootstrap) 
             # ignore original boostrap sample i and use instead new sample j:
             do.call(c, lapply(1:5, function(x) # resample only within quintiles, not across
               sample(all[d$quintile_Mbp5 == x], replace = T))))
  b = d[j, ] %>%
    group_by(., bin_Mbp5, quintile_Mbp5) %>%
    summarise(f4_num = sum(num_stat)/sum(num_sites),
              f4_denom = sum(denom_stat)/sum(denom_sites),
              f4_alpha = f4_num/f4_denom,
              f4_mex = 1 - f4_alpha) %>%
    arrange(quintile_Mbp5)
  return(c(spearman = with(b, cor(x = quintile_Mbp5, 
                                  y = f4_mex, 
                                  method = "spearman")),
           # also get alpha estimates for each quintile in the bootstrap sample
           unlist(pivot_wider(b, names_from = quintile_Mbp5, 
                              names_prefix = "Mbp", 
                              values_from = f4_mex, 
                              id_cols = quintile_Mbp5))))
}

# test:
#filter(d_f4, !is.na(num_sites)) %>%
#  calc_cor_r5(., 1:nrow(.))
#filter(d_f4, !is.na(num_sites) & !inv4m) %>%
#  calc_cor_r5(., 1:nrow(.))
#filter(d_f4, !is.na(num_sites)) %>%
#  calc_cor_r5(., 5)


# bootstrap resampling each recombination rate quintile independently
boot_r5 <- boot(data = filter(d_f4, !is.na(num_sites)), 
                statistic = calc_cor_r5,
                sim = "ordinary",
                stype = "i",
                R = n_boot)

#boot.ci(boot_r5, index = 1, type = c("norm", "basic", "perc"))

# bootstrap resampling each recombination rate quintile independently, BUT
# excluding all windows overlapping inv4m
boot_r5_no_inv4m <- boot(data = filter(d_f4, !is.na(num_sites) & !inv4m), 
                         statistic = calc_cor_r5,
                         sim = "ordinary",
                         stype = "i",
                         R = n_boot)

#boot.ci(boot_r5_no_inv4m, index = 1, type = c("norm", "basic", "perc"))



# bootstrap function to calc stat
# correlation with coding density (cd)
calc_cor_cd <- function(d, i) { # d is data, i is indices
  b = d[i, ] %>%
    group_by(., bin_cd5, quintile_cd5) %>%
    summarise(f4_num = sum(num_stat)/sum(num_sites),
              f4_denom = sum(denom_stat)/sum(denom_sites),
              f4_alpha = f4_num/f4_denom,
              f4_mex = 1 - f4_alpha) %>%
    arrange(quintile_cd5)
  return(c(spearman = with(b, cor(x = quintile_cd5, 
                                  y = f4_mex, 
                                  method = "spearman")),
           # also get estimates for each quintile in the bootstrap sample
           unlist(pivot_wider(b, names_from = quintile_cd5, 
                              names_prefix = "cd", 
                              values_from = f4_mex, 
                              id_cols = quintile_cd5))))
}


# test:
#filter(d_f4, !is.na(num_sites)) %>%
#  calc_cor_cd(., 1:nrow(.)) # all windows
#filter(d_f4, !is.na(num_sites) & !inv4m) %>%
#  calc_cor_cd(., 1:nrow(.)) # excluding inv4m

# bootstrap of correlation across coding density
boot_cd <- boot(data = filter(d_f4, !is.na(num_sites)), 
                statistic = calc_cor_cd,
                sim = "ordinary",
                stype = "i",
                R = n_boot)

# bootstrap ci for spearman's rank correlation:
#boot.ci(boot_cd, index = 1, type = c("norm", "basic", "perc"))

# bootstrap of correlation across coding density
boot_cd_no_inv4m <- boot(data = filter(d_f4, !is.na(num_sites) & !inv4m), 
                         statistic = calc_cor_cd,
                         sim = "ordinary",
                         stype = "i",
                         R = n_boot)

# bootstrap ci for spearman's rank correlation:
#boot.ci(boot_cd_no_inv4m, index = 1, type = c("norm", "basic", "perc"))



# make plots, add spearman's rank + basic bootstrap to both
# points and bars are from bootstap estimates
d_boot_r <- as.data.frame(rbind(boot_r$t0, boot_r$t)) %>%
  data.table::setnames(names(boot_r$t0)) %>%
  dplyr::mutate(boot = 0:n_boot) %>% # original estimate I set to boot '0'
  pivot_longer(., cols = starts_with("Mbp"), 
               names_to = "mbp_bin", values_to = "f4_mex")
ci_boot_r <- t(sapply(1:5, function(x)
  boot.ci(boot_r, index = 1+x, conf = 0.95, type = "basic")$basic[4:5])) %>%
  data.frame(.) %>%
  data.table::setnames(c("low", "high")) %>%
  mutate(mbp_bin = paste0("Mbp", 1:5))
ci_spearman_r = boot.ci(boot_r, index = 1, conf = 0.95, type = "basic")
text_spearman_r = paste0("Spearman's rho = ", ci_spearman_r$t0, 
                         " CI_95%(", ci_spearman_r$basic[4], ", ",
                         ci_spearman_r$basic[5], ")")

p_r5 <- ggplot(data = d_boot_r, aes(x = mbp_bin)) +
  # first plot original point estimates for ind. ancestry
  geom_point(data = filter(d_boot_r, boot > 0), 
             color = col_maize_mex_parv[zea],
             position = position_jitter(0.2),
             size = 0.1,
             aes(y = f4_mex)) +
  # then add mean for that group
  geom_point(data = filter(d_boot_r, boot == 0), 
             pch = 18, size = 2,
             aes(y = f4_mex)) +
  # and errorbars for 90% CI around that mean
  # based on bootstrap with NGSAdmix
  geom_errorbar(data = ci_boot_r, aes(ymin = low,
                                      ymax = high),
                width = .5) +
  xlab("Mbp per cM (quintiles low -> high)") +
  ylab("Proportion mexicana ancestry") +
  labs(subtitle = text_spearman_r) +
  theme_classic() +
  guides(color = guide_legend("Subspecies"),
         shape = guide_legend("Subspecies")) +
  ggtitle(paste("Ancestry by bp density in", sympatric_pop)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p_r5
ggsave(file = png_r5,
       plot = p_r5,
       device = "png",
       width = 5.5, height = 4.5, 
       units = "in", dpi = 300)

# plot w/out inv4m:
d_boot_r_no_inv4m <- as.data.frame(rbind(boot_r_no_inv4m$t0, boot_r_no_inv4m$t)) %>%
  data.table::setnames(names(boot_r_no_inv4m$t0)) %>%
  dplyr::mutate(boot = 0:n_boot) %>% # original estimate I set to boot '0'
  pivot_longer(., cols = starts_with("Mbp"), 
               names_to = "mbp_bin", values_to = "f4_mex")
ci_boot_r_no_inv4m <- t(sapply(1:5, function(x)
  boot.ci(boot_r_no_inv4m, index = 1+x, conf = 0.95, type = "basic")$basic[4:5])) %>%
  data.frame(.) %>%
  data.table::setnames(c("low", "high")) %>%
  mutate(mbp_bin = paste0("Mbp", 1:5))
ci_spearman_r_no_inv4m = boot.ci(boot_r_no_inv4m, index = 1, conf = 0.95, type = "basic")
text_spearman_r_no_inv4m = paste0("Spearman's rho = ", ci_spearman_r_no_inv4m$t0, 
                                  " CI_95%(", ci_spearman_r_no_inv4m$basic[4], ", ",
                                  ci_spearman_r_no_inv4m$basic[5], ")")

p_r5_no_inv4m <- ggplot(data = d_boot_r_no_inv4m, aes(x = mbp_bin)) +
  # first plot original point estimates for ind. ancestry
  geom_point(data = filter(d_boot_r_no_inv4m, boot > 0), 
             color = col_maize_mex_parv[zea],
             position = position_jitter(0.2),
             size = 0.1,
             aes(y = f4_mex)) +
  # then add mean for that group
  geom_point(data = filter(d_boot_r_no_inv4m, boot == 0), 
             pch = 18, size = 2,
             aes(y = f4_mex)) +
  # and errorbars for 90% CI around that mean
  # based on bootstrap with NGSAdmix
  geom_errorbar(data = ci_boot_r_no_inv4m, aes(ymin = low,
                                               ymax = high),
                width = .5) +
  xlab("Mbp per cM (quintiles low -> high)") +
  ylab("Proportion mexicana ancestry") +
  labs(subtitle = text_spearman_r_no_inv4m) +
  theme_classic() +
  guides(color = guide_legend("Subspecies"),
         shape = guide_legend("Subspecies")) +
  ggtitle(paste("Ancestry by bp density in", sympatric_pop, "no inv4m")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p_r5_no_inv4m
ggsave(file = png_r5_no_inv4m,
       plot = p_r5_no_inv4m,
       device = "png",
       width = 5.5, height = 4.5, 
       units = "in", dpi = 300)


# by coding bp per cM
d_boot_cd <- as.data.frame(rbind(boot_cd$t0, boot_cd$t)) %>%
  data.table::setnames(names(boot_cd$t0)) %>%
  dplyr::mutate(boot = 0:n_boot) %>% # original estimate I set to boot '0'
  pivot_longer(., cols = starts_with("cd"), 
               names_to = "cd_bin", values_to = "f4_mex")
ci_boot_cd <- t(sapply(1:5, function(x)
  boot.ci(boot_cd, index = 1+x, conf = 0.95, type = "basic")$basic[4:5])) %>%
  data.frame(.) %>%
  data.table::setnames(c("low", "high")) %>%
  mutate(cd_bin = paste0("cd", 1:5))
ci_spearman_cd = boot.ci(boot_cd, index = 1, conf = 0.95, type = "basic")
text_spearman_cd = paste0("Spearman's rho = ", ci_spearman_cd$t0, 
                          " CI_95%(", ci_spearman_cd$basic[4], ", ",
                          ci_spearman_cd$basic[5], ")")

p_cd5 <- ggplot(data = d_boot_cd, aes(x = cd_bin)) +
  # first plot original point estimates for ind. ancestry
  geom_point(data = filter(d_boot_cd, boot > 0), 
             color = col_maize_mex_parv[zea],
             position = position_jitter(0.2),
             size = 0.1,
             aes(y = f4_mex)) +
  # then add original estimate for that group
  geom_point(data = filter(d_boot_cd, boot == 0), 
             pch = 18, size = 2,
             aes(y = f4_mex)) +
  # and errorbars for 90% CI around that mean
  # based on bootstrap with NGSAdmix
  geom_errorbar(data = ci_boot_cd, aes(ymin = low,
                                       ymax = high),
                width = .5) +
  xlab("Coding bp per cM (quintiles low -> high)") +
  ylab("Proportion mexicana ancestry") +
  labs(subtitle = text_spearman_cd) +
  theme_classic() +
  guides(color = guide_legend("Subspecies"),
         shape = guide_legend("Subspecies")) +
  ggtitle(paste("Ancestry by gene density in", sympatric_pop)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p_cd5
ggsave(file = png_cd5,
       plot = p_cd5,
       device = "png",
       width = 5.5, height = 4.5, 
       units = "in", dpi = 300)


# by coding bp per cM, excluding windows overlapping inversion inv4m
d_boot_cd_no_inv4m <- as.data.frame(rbind(boot_cd_no_inv4m$t0, boot_cd_no_inv4m$t)) %>%
  data.table::setnames(names(boot_cd_no_inv4m$t0)) %>%
  dplyr::mutate(boot = 0:n_boot) %>% # original estimate I set to boot '0'
  pivot_longer(., cols = starts_with("cd"), 
               names_to = "cd_bin", values_to = "f4_mex")
ci_boot_cd_no_inv4m <- t(sapply(1:5, function(x)
  boot.ci(boot_cd_no_inv4m, index = 1+x, conf = 0.95, type = "basic")$basic[4:5])) %>%
  data.frame(.) %>%
  data.table::setnames(c("low", "high")) %>%
  mutate(cd_bin = paste0("cd", 1:5))
ci_spearman_cd_no_inv4m = boot.ci(boot_cd_no_inv4m, index = 1, conf = 0.95, type = "basic")
text_spearman_cd_no_inv4m = paste0("Spearman's rho = ", ci_spearman_cd_no_inv4m$t0, 
                                   " CI_95%(", ci_spearman_cd_no_inv4m$basic[4], ", ",
                                   ci_spearman_cd_no_inv4m$basic[5], ")")

p_cd5_no_inv4m <- ggplot(data = d_boot_cd_no_inv4m, aes(x = cd_bin)) +
  # first plot original point estimates for ind. ancestry
  geom_point(data = filter(d_boot_cd_no_inv4m, boot > 0), 
             color = col_maize_mex_parv[zea],
             position = position_jitter(0.2),
             size = 0.1,
             aes(y = f4_mex)) +
  # then add original estimate for that group
  geom_point(data = filter(d_boot_cd_no_inv4m, boot == 0), 
             pch = 18, size = 2,
             aes(y = f4_mex)) +
  # and errorbars for 90% CI around that mean
  # based on bootstrap with NGSAdmix
  geom_errorbar(data = ci_boot_cd_no_inv4m, aes(ymin = low,
                                                ymax = high),
                width = .5) +
  xlab("Coding bp per cM (quintiles low -> high)") +
  ylab("Proportion mexicana ancestry") +
  labs(subtitle = text_spearman_cd_no_inv4m) +
  theme_classic() +
  guides(color = guide_legend("Subspecies"),
         shape = guide_legend("Subspecies")) +
  ggtitle(paste("Ancestry by gene density in", sympatric_pop, "no inv4m")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p_cd5_no_inv4m
ggsave(file = png_cd5_no_inv4m,
       plot = p_cd5_no_inv4m,
       device = "png",
       width = 5.5, height = 4.5, 
       units = "in", dpi = 300)


# plot separately f4's for numerator and denominator
p_f4_num_denom <- bind_rows(f4_by_r, f4_by_cd) %>%
  pivot_longer(., cols = c("f4_num", "f4_denom"), 
               names_to = "num_denom", values_to = "f4") %>%
  pivot_longer(., cols = c("quintile_Mbp5", "quintile_cd5"),
               names_to = "type", values_to = "quintile") %>%
  filter(!is.na(quintile)) %>% # pivot longer creates some NAs because rows either have cd5 or Mbp5 data
  mutate(neg_f4 = -f4) %>%
  ggplot(., aes(x = quintile, y = neg_f4, color = num_denom)) +
  geom_point() +
  facet_wrap(~type) +
  ggtitle(paste("f4 estimates for", sympatric_pop, "(num) and allopatric maize (denom)")) +
  xlab("Mbp per cM quintile (low -> high)") +
  ylab("f4 estimate (negative of angsd output)") +
  labs(subtitle = "f4_num/f4_denom = alpha (maize ancestry estimate)",
       color = "f4 ratio components") +
  theme_light() +
  geom_hline(yintercept = 0)
# p_f4_num_denom
ggsave(file = png_f4_num_denom,
       plot = p_f4_num_denom,
       device = "png",
       width = 7.5, height = 4.5, 
       units = "in", dpi = 300)
library(dplyr)
library(tidyr)
library(ggplot2)
library(boot)
library(xtable)

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
zea = ifelse(sympatric_pop == "sympatric_maize" || (substr(sympatric_pop, 1, 3) == "pop" && as.integer(substr(sympatric_pop, 4, 7)) > 100), "maize", "mexicana") 

f4_num_file = snakemake@input[["f4_num"]]
# f4_num_file = paste0("ancestry_by_r/results/f4/", sympatric_pop, ".f4")
f4_denom_file = snakemake@input[["f4_denom"]]
# f4_denom_file = paste0("ancestry_by_r/results/f4/pop22.f4")
windows_file = snakemake@input[["windows"]]
# windows_file = "ancestry_by_r/results/map_pos_1cM_windows.txt"
colors_file = snakemake@input[["colors"]]
# colors_file = "colors.R"
png_r5 = snakemake@output[["png_r5"]]
# png_r5 = paste0("ancestry_by_r/plots/f4_pop22_symp_", sympatric_pop, "_byr5.png")
png_cd5 = snakemake@output[["png_cd5"]]
# png_cd5 = paste0("ancestry_by_r/plots/f4_pop22_symp_", sympatric_pop, "_bycd5.png")
png_r5_no_inv4m = snakemake@output[["png_r5_no_inv4m"]]
# png_r5_no_inv4m = paste0("ancestry_by_r/plots/f4_pop22_symp_", sympatric_pop, "_byr5_noinv4m.png")
png_cd5_no_inv4m = snakemake@output[["png_cd5_no_inv4m"]]
# png_cd5_no_inv4m = paste0("ancestry_by_r/plots/f4_pop22_symp_", sympatric_pop, "_bycd5_noinv4m.png")
png_f4_num_denom = snakemake@output[["png_f4_num_denom"]]
# png_f4_num_denom = paste0("ancestry_by_r/plots/f4_pop22_symp_", sympatric_pop, "_num_denom.png")
n_boot = snakemake@params[["n_boot"]]
# n_boot = 1000
inv_file = snakemake@input[["inv_file"]]
# inv_file = "data/refMaize/inversions/knownInv_v4_coord.txt"
file_pearsons_rho_f4 = snakemake@output[["tbl_pearsons_rho_f4"]]
file_pearsons_rho_f4 = paste0("ancestry_by_r/tables/tbl_pearsons_rho_f4_", zea, ".tex")

source(colors_file)

# f4-ratio calculation for alpha, admixture proportion from mexicana into X:
# f4(trip, parv; X, allo_maize) = E[(x_trip - x_parv)*(x_X - x_allo_maize)] in phylogeny (((parv, allo_maize), allo_mex), trip)
# where X is a sympatric population (numerator)
# or X is non-admixed allopatric mexicana (denominator)

# f4 ratio numerator
f4_num <- read.table(f4_num_file, stringsAsFactors = F, sep = "\t", header = T)
# f4 ratio denominator allopatric maize
f4_denom <- read.table(f4_denom_file, stringsAsFactors = F, sep = "\t", header = T)
# angsd outputs sum of (x_trip - x_parv)*(x_X - x_allo_maize)
# so I divide by number of sites to get the expectation
alpha = (sum(f4_num$Numer)/sum(f4_num$numSites))/(sum(f4_denom$Numer)/sum(f4_denom$numSites))
# genomewide the alpha estimate for % mexicana is ~ 21% in sympatric maize
# alpha

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


# combine data
d_f4 <- group_by(f4_num, window) %>%
  summarise(num_stat = sum(Numer),
            num_sites = sum(numSites)) %>%
  full_join(., group_by(f4_denom, window) %>%
              summarise(denom_stat = sum(Numer),
                        denom_sites = sum(numSites)),
            by = "window") %>%
  left_join(winds, ., by = "window") %>%
  ungroup()

# by recombination rate
f4_by_r <- d_f4 %>%
  group_by(., bin_r5, quintile_r5) %>%
  filter(!is.na(num_sites)) %>% # eliminates 2 windows w/out data
  summarise(f4_num = sum(num_stat)/sum(num_sites),
            f4_denom = sum(denom_stat)/sum(denom_sites),
            f4_alpha = f4_num/f4_denom,
            n_sites_tot = (sum(num_sites) + sum(denom_sites))/2,
            n_windows = length(unique(window))) %>%
  arrange(quintile_r5) %>%
  ungroup()
#View(f4_by_r)

# by coding density
f4_by_cd <- d_f4 %>%
  group_by(., bin_cd5, quintile_cd5) %>%
  filter(!is.na(num_sites)) %>% # eliminates 2 windows w/out data
  summarise(f4_num = sum(num_stat)/sum(num_sites),
            f4_denom = sum(denom_stat)/sum(denom_sites),
            f4_alpha = f4_num/f4_denom,
            n_sites_tot = (sum(num_sites) + sum(denom_sites))/2,
            n_windows = length(unique(window))) %>%
  arrange(quintile_cd5) %>%
  ungroup()
#View(f4_by_cd)


# bootstrap

# re-sampling is done within each quintile separately
# bootstrapping by recombination rate function
calc_cor_r5 <- function(d, i) { # d is data, i is indices
  all = 1:nrow(d)
  j = ifelse(i == all, i, # original sample (not bootstrap) 
             # ignore original boostrap sample i and use instead new sample j:
             do.call(c, lapply(1:5, function(x) # resample only within quintiles, not across
               sample(all[d$quintile_r5 == x], replace = T))))
  b = d[j, ] %>%
    group_by(., bin_r5, quintile_r5) %>%
    summarise(f4_num = sum(num_stat)/sum(num_sites),
              f4_denom = sum(denom_stat)/sum(denom_sites),
              f4_alpha = f4_num/f4_denom) %>%
    arrange(quintile_r5)
  return(c(spearman = with(b, cor(x = quintile_r5, 
                                  y = f4_alpha, 
                                  method = "spearman")),
           # also get alpha estimates for each quintile in the bootstrap sample
           unlist(pivot_wider(b, names_from = quintile_r5, 
                              names_prefix = "r", 
                              values_from = f4_alpha, 
                              id_cols = quintile_r5))))
}

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


# re-sampling is done within each quintile separately
# bootstrapping by coding density (cd) function:
calc_cor_cd5 <- function(d, i) { # d is data, i is indices
  all = 1:nrow(d)
  j = ifelse(i == all, i, # original sample (not bootstrap) 
             # ignore original boostrap sample i and use instead new sample j:
             do.call(c, lapply(1:5, function(x) # resample only within quintiles, not across
               sample(all[d$quintile_cd5 == x], replace = T))))
  b = d[j, ] %>%
    group_by(., bin_cd5, quintile_cd5) %>%
    summarise(f4_num = sum(num_stat)/sum(num_sites),
              f4_denom = sum(denom_stat)/sum(denom_sites),
              f4_alpha = f4_num/f4_denom) %>%
    arrange(quintile_cd5)
  return(c(spearman = with(b, cor(x = quintile_cd5, 
                                  y = f4_alpha, 
                                  method = "spearman")),
           # also get estimates for each quintile in the bootstrap sample
           unlist(pivot_wider(b, names_from = quintile_cd5, 
                              names_prefix = "cd", 
                              values_from = f4_alpha, 
                              id_cols = quintile_cd5))))
}

# bootstrap of correlation across coding density
boot_cd5 <- boot(data = filter(d_f4, !is.na(num_sites)), 
                statistic = calc_cor_cd5,
                sim = "ordinary",
                stype = "i",
                R = n_boot)

# bootstrap ci for spearman's rank correlation:
#boot.ci(boot_cd5, index = 1, type = c("norm", "basic", "perc"))

# bootstrap of correlation across coding density
boot_cd5_no_inv4m <- boot(data = filter(d_f4, !is.na(num_sites) & !inv4m), 
                         statistic = calc_cor_cd5,
                         sim = "ordinary",
                         stype = "i",
                         R = n_boot)

# bootstrap ci for spearman's rank correlation:
#boot.ci(boot_cd5_no_inv4m, index = 1, type = c("norm", "basic", "perc"))

# make 4 plots for 4 different bootstrap analyses (cd/r * with/without inv4m):
plots = c(png_r5, png_r5_no_inv4m, png_cd5, png_cd5_no_inv4m)
boots = list(boot_r5, boot_r5_no_inv4m, boot_cd5, boot_cd5_no_inv4m)
feature = c("r", "r", "cd", "cd")
x_axis_labels_r = filter(winds, !duplicated(paste(bin_r5))) %>%
  arrange(quintile_r5) %>%
  dplyr::select(bin_r5) %>%
  rename(bin = bin_r5)
x_axis_labels_cd = filter(winds, !duplicated(paste(bin_cd5))) %>%
  arrange(quintile_cd5) %>%
  dplyr::select(bin_cd5) %>%
  rename(bin = bin_cd5)
x_axis_labels = c(x_axis_labels_r, x_axis_labels_r, x_axis_labels_cd, x_axis_labels_cd)
x_labs = c("Recombination rate quintile (cM/Mb)", "Recombination rate quintile (cM/Mb)",
           "Gene density quintile (coding bp/cM)", "Gene density quintile (coding bp/cM)")
for (i in 1:4){
  # make plots
  d_boot <- as.data.frame(rbind(boots[[i]]$t0, boots[[i]]$t)) %>%
    data.table::setnames(names(boots[[i]]$t0)) %>%
    dplyr::mutate(boot = 0:n_boot) %>% # original estimate I set to boot '0'
    pivot_longer(., cols = starts_with(feature[i]), 
                 names_to = "bin", 
                 values_to = "f4_alpha")
  ci_boot <- t(sapply(1:5, function(x)
    boot.ci(boots[[i]], index = 1 + x, conf = 0.95, type = "perc")$perc[4:5])) %>%
    data.frame(.) %>%
    data.table::setnames(c("low", "high")) %>%
    mutate(bin = paste0(feature[i], 1:5))
  
  p <- ggplot(data = d_boot, aes(x = bin)) +
    # violin plot of individual bootstrap samples
    geom_violin(data = filter(d_boot, boot > 0),
                color = col_maize_mex_parv[zea],
                fill = col_maize_mex_parv[zea],
                aes(y = f4_alpha)) +
    # then add mean from original estimate (not bootstraps) 
    geom_point(data = filter(d_boot, boot == 0),
               pch = 18, size = 2,
               aes(y = f4_alpha)) +
    # and errorbars for 95% CI around that mean
    geom_errorbar(data = ci_boot, aes(ymin = low,
                                      ymax = high),
                  width = .2) +
    xlab(x_labs[i]) +
    ylab("Proportion mexicana ancestry") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = x_axis_labels[i]$bin) +
    ylim(0:1)
  #p
  ggsave(file = plots[i],
         plot = p,
         device = "png",
         width = 5.5, height = 4.5,
         units = "in", dpi = 300)
}

# plot separately f4's for numerator and denominator
p_f4_num_denom <- bind_rows(f4_by_r, f4_by_cd) %>%
  pivot_longer(., cols = c("f4_num", "f4_denom"),
               names_to = "num_denom", values_to = "f4") %>%
  pivot_longer(., cols = c("quintile_r5", "quintile_cd5"),
               names_to = "type", values_to = "quintile") %>%
  filter(!is.na(quintile)) %>% # pivot longer creates some NAs because rows either have cd5 or Mbp5 data
  mutate(neg_f4 = -f4) %>%
  ggplot(., aes(x = quintile, y = neg_f4, color = num_denom)) +
  geom_point() +
  facet_wrap(~type) +
  ggtitle(paste("f4 estimates for sympatric", zea, "(num) and allopatric maize (denom)")) +
  xlab("Recombination rate quintile (cM/Mb)") +
  ylab("f4 estimate (negative of angsd output)") +
  labs(subtitle = "alpha = f4_num/f4_denom",
       color = "f4 ratio components") +
  theme_light() +
  geom_hline(yintercept = 0) +
  scale_x_discrete(labels = x_axis_labels_r$bin_r5)
# p_f4_num_denom
ggsave(file = png_f4_num_denom,
       plot = p_f4_num_denom,
       device = "png",
       width = 7.5, height = 4.5,
       units = "in", dpi = 300)

# make table summarising pearson's correlations:
rho = data.frame(method = "f4 ratio",
                 feature = c(rep("recombination rate (cM/Mb)", 2),
                             rep("gene density (coding bp/cM)", 2)),
                 resolution = "genomic quintiles",
                 data = rep(c("whole genome", "excl. inv4m"), 2),
                 group = paste("sympatric", zea),
                 stringsAsFactors = F) %>%
  mutate(
    rho_estimate = sapply(1:4, function(i) boots[[i]]$t0[["spearman"]]),
    boot_low = sapply(1:4, function(i) boot.ci(boots[[i]], index = 1, 
                                               conf = 0.95, type = "perc")$perc[4]),
    boot_high = sapply(1:4, function(i) boot.ci(boots[[i]], index = 1, 
                                               conf = 0.95, type = "perc")$perc[5])) %>%
  #dplyr::select(method, feature, resolution, group, rho_estimate, boot_low, boot_high) %>%
  rename(`Pearson's rank correlation` = rho_estimate, `2.5%` = boot_low, `97.5%` = boot_high)

# print table to file for estimates of Pearson's rank correlation
print(xtable(rho, 
             digits = 3,
             label = paste0("tbl_pearsons_rho_f4_", zea),
             type = "latex", 
             latex.environments = NULL),
      include.rownames = F,
      file = file_pearsons_rho_f4)



# This script plots the the bootstrap results for 
# ancestry across recombination quintiles
library(dplyr)
library(tidyr)
library(ggplot2)

PREFIX <- "pass2_alloMAIZE_PalmarChico"
WIND = "1cM"
K = 3
PREFIX_METRICS <- "hilo_alloMAIZE_MAIZE4LOW"
ancestries <- c("maize", "mexicana", "parviglumis")

# get meta data
metrics <- read.table(paste0("../filtered_bams/metrics/", PREFIX_METRICS, ".flagstat.total"), header = F, stringsAsFactors = F)
colnames(metrics) <- c("ID", "total_reads_pass")
hilo <- read.table("../samples/hilo_meta.txt", stringsAsFactors = F, header = T, sep = "\t")
landraces <- read.table("../samples/alloMAIZE_meta.txt", stringsAsFactors = F, header = T, sep = "\t")
parviglumis <- read.table("../samples/PalmarChico_meta.txt", stringsAsFactors = F, header = T, sep = "\t")

# a rough coverage estimate is number of reads passing filtering * 150bp/read
# divided by total area they could map to in maize reference genome v4

# size of reference genome reads are mapped to
ref_genome_size <- sum(read.table("../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.fai", 
                                  stringsAsFactors = F)$V2) # V2 is the size of each chromosome mapped to,
# including parts I won't analyze on the Pt and Mt and scaffolds not assigned to chromosomes (but reads would pass Q filters there too)
metrics$est_coverage = round(metrics$total_reads_pass*150/ref_genome_size, 4)

# combine sample meta data
meta <- bind_rows(hilo, landraces, parviglumis) %>%
  left_join(., metrics, by = "ID") %>%
  mutate(., group = paste(symp_allo, zea, sep = "_"))

# get bootstrap data of ancestry %
anc_boot <- do.call(rbind,
               lapply(0:100, function(BOOT) do.call(rbind,
                                          lapply(1:5, function(r) 
                                            read.table(paste0("results/bootstrap/windows_", WIND,
                         "/r5_recomb", r, "/", PREFIX, "/K", K, "/boot", BOOT, ".anc"),
                  header = T, sep = "\t") %>%
  mutate(bootstrap = BOOT) %>%
  mutate(r5_bin = r))))) %>%
  left_join(., meta, by = "ID")

# set sequencing coverage cutoff for inclusion
min_coverage = 0.1

anc_boot_mean <- anc_boot %>%
  filter(zea == "parviglumis" | 
           est_coverage >= min_coverage) %>%
  group_by(group, zea, symp_allo, bootstrap, r5_bin) %>%
  summarise(mexicana_ancestry = mean(mexicana),
            maize_ancestry = mean(maize),
            parviglumis_ancestry = mean(parviglumis))
# point estimate individual ancestry
anc_ind_mean <- anc_boot %>%
  filter(zea == "parviglumis" | 
           est_coverage >= min_coverage) %>%
  filter(bootstrap == 1) %>% 
  gather(., key = "ancestry", value = "p", ancestries)
# bootstrap 90% confidence interval individual ancestry
anc_ind_boot_confidence_intervals <- anc_boot %>%
  filter(bootstrap != 0) %>% # this is the original sample
  gather(., key = "ancestry", value = "p", 
         ancestries) %>%
  group_by(ID, ancestry, r5_bin) %>%
  summarise(low_boot = nth(p, 6, order_by = p), # low and high bounds
            high_boot = nth(p, 95, order_by = p), # of 90% conf. interval
            mid_boot = nth(p, 50, order_by = p),
            mean_boot = mean(p))
anc_ind_estimate <- anc_ind_mean %>%
  left_join(., anc_ind_boot_confidence_intervals, by = c("ancestry", "ID", "r5_bin")) %>% 
  mutate(low = low_boot - mean_boot + p,
         high = high_boot - mean_boot + p)


# plot
anc_boot_mean %>%
  ggplot(., aes(group = r5_bin, y = maize_ancestry)) +
  geom_boxplot() +
  facet_wrap(~group) +
  ggtitle("Mean maize ancestry K3 bootstrapped")
anc_boot_mean %>%
  ggplot(., aes(group = r5_bin, y = mexicana_ancestry)) +
  geom_boxplot() +
  facet_wrap(~group) +
  ggtitle("Mean mexicana ancestry K3 bootstrapped")
anc_boot_mean %>%
  ggplot(., aes(group = r5_bin, y = parviglumis_ancestry)) +
  geom_boxplot() +
  facet_wrap(~group) +
  ggtitle("Mean parviglumis ancestry K3 bootstrapped")

# plot all together
anc_boot_mean %>%
  gather(., key = "ancestry", value = "p", 
         paste(ancestries, "ancestry", sep = "_")) %>%
  ggplot(aes(x = as.factor(r5_bin), y = p, color = ancestry)) +
  geom_boxplot() +
  facet_wrap(~group) +
  xlab("recombination quintile (low to high r)") +
  ggtitle("Mean K3 ancestries bootstrapped")

anc_boot_mean %>%
  ggplot(aes(x = as.factor(r5_bin), y = mexicana_ancestry,
             color = group)) +
  geom_boxplot(outlier.size = 1, notch = T) +
  geom_jitter(shape = 16, size = .1, position=position_jitter(0.2)) +
  #geom_point(size = .1, color = "black") +
  facet_wrap(~group) +
  ggtitle("Mean K3 mexicana ancestry bootstrapped")

anc_group_mean <- anc_boot_mean %>%
  filter(bootstrap == 0) %>%
  gather(., key = "ancestry", value = "p", 
         paste(ancestries, "ancestry", sep = "_"))
anc_boot_confidence_intervals <- anc_boot_mean %>%
  filter(bootstrap != 0) %>% # this is the original sample
  gather(., key = "ancestry", value = "p", 
         paste(ancestries, "ancestry", sep = "_")) %>%
  group_by(ancestry, group, r5_bin) %>%
  summarise(low_boot = nth(p, 6, order_by = p), # low and high bounds
            high_boot = nth(p, 95, order_by = p), # of 90% conf. interval
            mid_boot = nth(p, 50, order_by = p),
            mean_boot = mean(p))
anc_group_estimate <- anc_group_mean %>%
  left_join(., anc_boot_confidence_intervals, by = c("ancestry", "group", "r5_bin")) %>%
  mutate(low = low_boot - mean_boot + p,
         high = high_boot - mean_boot + p)

# plot 90% confidence interval around mean from bootstrap:
anc_group_estimate %>%
  ggplot(aes(x = as.factor(r5_bin), y = p, 
             color = ancestry)) + 
  geom_errorbar(aes(ymin = low, 
                    ymax = high), 
                width = .5) +
  geom_point(pch = 18, size = 1) +
  facet_wrap(~group) +
  xlab("recombination quintile (low to high r)") +
  ggtitle("Mean K3 ancestry estimates and 90% bootstrap CI")
ggsave("plots/K3_mean_ancestry_all_bootstrap_90CI.png",
       width = 8, height = 8, units = "in")

# look at just mexicana ancestry in sympatric maize and mexicana:
anc_group_estimate %>%
  filter(zea != "parviglumis") %>%
  filter(ancestry == "mexicana_ancestry") %>%
  ggplot(aes(x = as.factor(r5_bin), y = p, group = zea)) +
  # first plot original point estimates for ind. ancestry
  geom_jitter(data = filter(anc_ind_mean, zea != "parviglumis" & ancestry == "mexicana") %>%
                mutate(depth = ifelse(est_coverage > 2, 2, est_coverage)),
              aes(x = as.factor(r5_bin),
                  y = p,
                  color = zea,
                  size = depth),
              shape = 16, 
              position = position_jitter(0.2)) +
  scale_size(range = c(.05, 2)) +
  # then add mean for that group
  geom_point(pch = 18, size = 2) +
  # and errorbars for 90% CI around that mean 
  # based on bootstrap with NGSAdmix
  geom_errorbar(aes(ymin = low, 
                    ymax = high), 
                width = .5) +
  facet_wrap(~symp_allo) +
  xlab("recombination quintile (low to high r)") +
  ggtitle("Ind. mexicana ancestry estimates (with mean and 90% bootstrap CI)")
ggsave("plots/K3_ind_ancestries_jitter_and_mean_bootstrap_90CI.png",
       width = 8, height = 8, units = "in")

# add jitter to individual ancestries & CI high and low bounds
anc_ind_jitter <- anc_ind_estimate %>%
  filter(., zea != "parviglumis" & ancestry == "mexicana") %>%
  mutate(r5_bin_jitter = jitter(r5_bin, factor = 2))
# plot with individual ancestry error bars
anc_group_estimate %>%
  filter(zea != "parviglumis") %>%
  filter(ancestry == "mexicana_ancestry") %>%
  ggplot(aes(x = r5_bin, y = p, group = zea)) +
  # first plot original point estimates for ind. ancestry
  geom_point(data = anc_ind_jitter,
                aes(x = r5_bin_jitter,
                    color = zea,
                    y = p)) +
  # add individual means error bars
  geom_errorbar(data = anc_ind_jitter,
                aes(x = r5_bin_jitter,
                    color = zea,
                    ymin = low,
                    ymax = high),
                width = .2) +
  # then add mean for that group
  geom_point(pch = 18, size = 2) +
  # and errorbars for 90% CI around that mean 
  # based on bootstrap with NGSAdmix
  geom_errorbar(aes(ymin = low, 
                    ymax = high), 
                width = .5) +
  facet_wrap(~symp_allo) +
  scale_y_continuous(breaks = seq(0, 1.2, by = .2)) +
  xlab("recombination quintile (low to high r)") +
  ggtitle("Ind. mexicana ancestry estimates (with mean and 90% bootstrap CI)")
ggsave("plots/K3_ind_ancestries_90CI_and_jitter_and_mean_bootstrap_90CI.png",
       width = 16, height = 8, units = "in")

# plot ancestry across recombination bins and elevation:
anc_ind_estimate %>%
  left_join(., dplyr::select(elev, popN, ELEVATION), by = "popN") %>%
  filter(., symp_allo == "sympatric") %>%
  filter(., ancestry == "mexicana") %>%
  ggplot(., aes(x = ELEVATION, y = p, 
                color = LOCALITY, 
                size = est_coverage, 
                shape = zea)) +
  geom_point() +
  ylab("mexicana ancestry") +
  geom_smooth(method = "lm", aes(group = zea), color = "black") +
  ggtitle("Higher mexicana ancestry at higher elevations") +
  facet_wrap(~r5_bin)
ggsave(paste0("plots/lm_predict_NGSadmix_proportion_mexicana-like_by_elevation_andByRecombBin_colorByPop_K", K, "_", PREFIX, ".png"), 
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)
ggsave(paste0("../../hilo_manuscript/figures/lm_predict_NGSadmix_proportion_mexicana-like_by_elevation_andByRecombBin_colorByPop_K", K, "_", PREFIX, ".png"), 
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)





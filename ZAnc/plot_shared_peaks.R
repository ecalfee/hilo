#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)

# this script plots ancestry outliers across the genome,
# from individual populations, and shared across pops

# load variables from Snakefile
zea = snakemake@params[["zea"]]
# zea = "maize"
colors_file = snakemake@input[["colors"]]
# colors_file = "colors.R"
sim_file = snakemake@input[["sim"]]
# sim_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".MVN.RData")
anc_file = snakemake@input[["anc"]]
# anc_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pops.anc.RData")
meta_file = snakemake@input[["meta_pop"]]
# meta_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pop.meta.RData")
sites_file = snakemake@input[["sites"]]
# sites_file = "local_ancestry/results/thinnedSNPs/HILO_MAIZE55/whole_genome.var.sites"
inv_file = snakemake@input[["inv"]]
# inv_file = "data/refMaize/inversions/knownInv_v4_coord.txt"
Ne = as.integer(snakemake@params[["Ne"]])
# Ne = 10000
YESNO = snakemake@params[["YESNO"]]
# YESNO = "yes"
prefix_all = snakemake@params[["prefix_all"]]
# prefix_all = "HILO_MAIZE55"

png_hist = paste0("ZAnc/plots/Ne", Ne, "_", YESNO, "Boot/", zea, "_hist_outlier_peaks.png")
png_ratio = paste0("ZAnc/plots/Ne", Ne, "_", YESNO, "Boot/", zea, "_ratio_outlier_peaks.png")
png_hist_ratio = paste0("ZAnc/plots/Ne", Ne, "_", YESNO, "Boot/", zea, "_hist_ratio_outlier_peaks.png")
png_chr_i_prefix = paste0("ZAnc/plots/Ne", Ne, "_", YESNO, "Boot/", zea, "_shared_outliers_chr_")

# load data
source(colors_file)
load(anc_file)
load(sim_file)
load(meta_file)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(cbPalette) <- c("black", "maize", "cb3", "cb4", "cb5", "cb6", "cb7", "mexicana")
if (zea == "maize"){
  col_obs_sim = c(cbPalette[zea], col_maize_mex_parv[[zea]])
}else{
  col_obs_sim = c(col_maize_mex_parv[[zea]], "plum")
}
names(col_obs_sim) = c("data", "mvn_sim")

# load inversion coordinates
inv = read.table(inv_file, stringsAsFactors = F, header = F) %>%
  data.table::setnames(c("ID", "chr", "start", "end", "length"))


# and site/position information for SNPs with ancestry calls
sites <- read.table(sites_file, header = F, stringsAsFactors = F,
                    sep = "\t") %>%
  data.table::setnames(c("chr", "pos", "major", "minor")) %>%
  dplyr::mutate(inv4m = (chr == inv$chr[inv$ID=="inv4m"] & 
             pos >= inv$start[inv$ID=="inv4m"] & 
             pos <= inv$end[inv$ID == "inv4m"])) # mark sites within inversion

# take inverse if mexicana to get 'introgressed' minor ancestry
if (zea == "mexicana"){
  anc = 1 - anc
  MVN_sim = 1 - MVN_sim
  meta_pops$alpha_local_ancestry = 1 - meta_pops$alpha_local_ancestry
}

# what is 2sd above the mean introgressed ancestry threshold for each pop?
meta_pops$sd2 <- apply(anc, 2, mean) + 2*apply(anc, 2, sd)

# find pop outliers in observed data
anc_outliers <- data.frame(anc, stringsAsFactors = F) %>%
  cbind(., sites) %>%
  tidyr::pivot_longer(., cols = colnames(anc), names_to = "pop", values_to = "anc") %>%
  left_join(., meta_pops, by = "pop") %>%
  mutate(top_sd2 = anc > sd2) %>%
  arrange(ELEVATION)

# find pop outliers in simulated data
sim_outliers <- data.frame(MVN_sim, stringsAsFactors = F) %>%
  mutate(chr = NA, # add placeholder position data
         pos = 1:nrow(.)) %>%
  tidyr::pivot_longer(., cols = colnames(anc), names_to = "pop", values_to = "anc") %>%
  left_join(., meta_pops, by = "pop") %>%
  mutate(top_sd2 = anc > sd2) %>%
  arrange(ELEVATION)

outliers_bypop <- bind_rows(mutate(anc_outliers, source = "data"),
                            mutate(sim_outliers, source = "mvn_sim", inv4m = F)) %>%
  group_by(chr, pos, source, inv4m) %>%
  summarise(outlier_pops = sum(top_sd2)) %>%
  filter(outlier_pops > 0) # only include if at least one outlier

# get ratio observed/simulated shared outliers:
# but first find a maximum # sharing pops with simulated values
# because otherwise the ratio is dividing by zero
# the max bin is 1 minus the first number of pops with no simulated shared outliers
max_bin = which.min(1:14 %in% unique(outliers_bypop$outlier_pops[outliers_bypop$source == "mvn_sim"])) - 1
# calculate ratio obs/sim outliers for every bin
outliers_ratios = outliers_bypop %>%
  dplyr::mutate(outlier_pops = ifelse(outlier_pops >= max_bin, max_bin, outlier_pops)) %>% # group outliers shared by >= max_bin pops so all bins have at least 1 count
  dplyr::group_by(source, outlier_pops) %>%
  dplyr::summarise(n = n()) %>% # count number of outliers each type (i.e. count of # sharing pops)
  tidyr::pivot_wider(data = ., 
                     names_from = source,
                     values_from = n) %>%
  dplyr::mutate(mvn_sim = ifelse(is.na(mvn_sim), 0, mvn_sim)) %>% # if no outliers found for a type, set to 0 (not NA)
  dplyr::mutate(data = data/sum(data), # get proportions from counts
                mvn_sim = mvn_sim/sum(mvn_sim),
                ratio = data/mvn_sim)
  

# plot sfs of shared outliers across pops
p_hist = outliers_bypop %>%
  ggplot(., aes(x = outlier_pops, stat(density))) +
  geom_histogram(aes(fill = source), 
                 position = "dodge2", alpha = .5, binwidth = 0.5) + 
  geom_histogram(data = filter(outliers_bypop, !inv4m),
                 aes(fill = source),
                 position = "dodge2",
                 binwidth = 0.5) +
  scale_fill_manual(values = col_obs_sim,
                    name = NULL, 
                    labels = c("observed", "simulated")) +
  ylim(c(0, 1.25)) +
  theme_classic() +
  xlab("Number of populations with high introgression") +
  ylab("Density") +
  ggtitle(paste("Distribution of shared outliers in", zea))
#p_hist
ggsave(file = png_hist,
       plot = p_hist,
       device = "png",
       height = 4, width = 5, 
       units = "in", dpi = 300)

# plot distribution of shared outliers a different way
# rather than histogram, show the ratio observed/simulated

p_ratio <- outliers_ratios %>% # get ratio observed/expected
  ggplot(., aes(x = outlier_pops, y = ratio, size = data)) +
  geom_point(color = col_maize_mex_parv[[zea]]) + 
  geom_line(size = 1, color = col_maize_mex_parv[[zea]]) + 
  theme_classic() +
  geom_abline(slope = 0, intercept = 1) +
  #ylim(c(0, 10)) +
  ylab("Ratio observed/expected outlier loci") +
  xlab("Number of populations with high introgression") +
  ggtitle(paste("Distribution of shared outliers across sympatric", zea, "populations")) +
  labs(size = "Proportion of \ntotal outliers") +
  scale_x_continuous(breaks = 1:max_bin, 
                     labels = c(1:(max_bin - 1), 
                                paste0(">=", max_bin)))
#p_ratio
ggsave(file = png_ratio,
       plot = p_ratio,
       height = 4, width = 6, units = "in",
       device = "png") 
  
# plot histogram and ratio on same plot with 2 axes:
scale_y2 = 5
p_hist_ratio = p_hist + 
  geom_point(data = outliers_ratios, aes(x = outlier_pops, y = ratio/scale_y2)) + # plot ratio observed/expected
  geom_line(data = outliers_ratios, aes(x = outlier_pops, y = ratio/scale_y2)) +
  geom_abline(slope = 0, intercept = 1/scale_y2, lty = "dashed") +
  scale_y_continuous(
    limits = c(0, 1.25),
    # Add a second axis based on a linear transformation of the first y axis
    sec.axis = sec_axis(trans=~.*scale_y2, name="Ratio observed/expected")
  ) +
  scale_fill_manual(values = col_obs_sim,
                    name = NULL, 
                    labels = c("observed", "expected"))
# p_hist_ratio
ggsave(file = png_hist_ratio,
       plot = p_hist_ratio,
       height = 4, width = 5, units = "in",
       device = "png") 

# make a plot of outliers on each chromosome
for (i in 1:10){
  p_chr_i = anc_outliers %>%
    filter(chr == i) %>%
    ggplot(., aes(pos/10^6, anc, color = top_sd2)) + # plot by position on chromosome (Mb), not relative position
    geom_point(size = .1) +
    geom_hline(data = meta_pops, 
               aes(yintercept = alpha_local_ancestry), 
               linetype = "dashed", alpha = .5,
               color = col_maize_mex_parv[[zea]]) +
    facet_grid(reorder(LOCALITY, desc(ELEVATION)) ~ .) +
    xlab("position (Mbp)") +
    ylab("introgressed ancestry frequency") +
    ggtitle(paste("Chr", i)) +
    ylim(c(0,1)) +
    scale_colour_manual(values = c("grey", "blue")) + 
    theme_classic() +
    coord_cartesian(ylim=c(0,1)) + 
    theme(legend.position = "none") +
    theme(strip.text.y = element_text(angle=0))
  ggsave(file = paste0(png_chr_i_prefix, i, ".png"),
         plot = p_chr_i,
         height = 8, width = 16, units = "in",
         device = "png")
}

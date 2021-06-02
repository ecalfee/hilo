#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
# this script plots timing of admixture estimated from ancestry_hmm

# load variables from Snakefile
Ne = snakemake@params[["Ne"]]
# Ne = 10000
YESNO = snakemake@params[["YESNO"]]
# YESNO = "yes"
alpha = as.numeric(snakemake@params[["alpha"]])
# alpha = 0.05
prefix = snakemake@params[["prefix_all"]]
# prefix = "HILO_MAIZE55_PARV50"
png_times = snakemake@output[["png_times"]]
# png_times = paste0("local_ancestry/plots/", prefix, "K3_admix_times_Ne", Ne, "_", YESNO, "Boot.png")
png_times_lzw = snakemake@output[["png_times_lzw"]]
# png_times_lzw = paste0("../hilo_manuscript/figures_supp/", prefix, "K3_admix_times_Ne", Ne, "_", YESNO, "Boot.tif")
txt_times = snakemake@output[["txt_times"]]
# txt_times = paste0("local_ancestry/results/", prefix, "K3_admix_times_Ne", Ne, "_", YESNO, "Boot.txt")
rds_times = snakemake@output[["rds_times"]]
# rds_times = paste0("local_ancestry/results/", prefix, "K3_admix_times_Ne", Ne, "_", YESNO, "Boot.RDS")
source(snakemake@params[["colors"]]) # plotting colors
# source("colors.R")

# load meta and bootstrap data for maize and mexicana
zea = c("maize", "mexicana")
for (z in zea){
  parv_maize = c("parv", "maize")
  load(paste0("local_ancestry/results/ancestry_hmm/",
              prefix, "/K3/Ne", Ne, "_", YESNO,
              "Boot/anc/", z, ".pop.meta.RData"))
  times = lapply(parv_maize, function(a)
                 do.call(bind_cols, lapply(meta_pops$pop, function(p)
    read.table(paste0("local_ancestry/results/ancestry_hmm/",
                      prefix, "/K3/Ne", Ne, "_", YESNO,
                      "Boot/", p, ".times"),
               header = T, stringsAsFactors = F, sep = " ") %>%
      dplyr::filter(ancestry == a) %>%
      dplyr::select(t))) %>%
      data.table::setnames(meta_pops$pop))
  names(times) = parv_maize
  boots = lapply(times, function(t) t[2:nrow(t), ]) # first row is estimate, remaining rows bootstraps
  names(boots) = parv_maize
  time_est = do.call(rbind,
                     lapply(parv_maize, function(a)
    data.frame(time = unlist(times[[a]][1, ]),
                    low_boot = apply(boots[[a]], 2, function(x) quantile(x, alpha/2)),
                    high_boot = apply(boots[[a]], 2, function(x) quantile(x, 1 - alpha/2)),
                    mean_boot = apply(boots[[a]], 2, mean),
                    sd_boot = apply(boots[[a]], 2, sd),
                    pop = colnames(times[[a]]), stringsAsFactors = F,
                    alpha = alpha) %>%
    dplyr::mutate(low_basic = 2*time - high_boot,
                  high_basic = 2*time - low_boot,
                  se_boot = sd_boot/sqrt(nrow(boots[[a]])),
                  admixture_pulse = a)))
  assign(x = z, value = left_join(meta_pops, time_est, by = "pop"))
  rm(meta_pops, times, time_est, boots)
}

# combine data
d <- bind_rows(maize, mexicana) %>%
  dplyr::mutate(introgress = ifelse(admixture_pulse == "maize", alpha_local_ancestry_maize, alpha_local_ancestry_parv))
  #dplyr::mutate(introgress = ifelse(zea == "mexicana", 1 - alpha_local_ancestry_mexicana, 1 - alpha_local_ancestry_maize))

# plot
p_times <- d %>%
  mutate(zea = paste("sympatric", zea)) %>%
  mutate(admixture_pulse = ifelse(admixture_pulse == "parv", "parviglumis", admixture_pulse)) %>%
  ggplot(., aes(x = reorder(LOCALITY, ELEVATION),
                y = time,
                shape = admixture_pulse,
                alpha = introgress < .1,
                #alpha = (introgress < .1 | introgress > .9),
                col = admixture_pulse)) +
  geom_point() + # plot time estimate
  # add errorbars for 95% percentile CI around that mean
  # based on bootstrap
  geom_errorbar(aes(ymin = low_boot,
                    ymax = high_boot),
                width = .5) +
  xlab("Sympatric population") +
  ylab("Admixture time (generations)") +
  theme_classic() +
  ylim(-100, 1650) +
  guides(color = guide_legend("Admixture pulse"),
         shape = guide_legend("Admixture pulse"),
         alpha = guide_legend("Pulse size")) +
  scale_color_manual(values = col_maize_mex_parv) +
  scale_shape_manual(values = c(1, 0)) +
  scale_alpha_manual(values = c(1, .3), labels = c("> 10% total ancestry", "< 10% total ancestry")) +
  #scale_alpha_manual(values = c(1, .3), labels = c("Intermediate\n(10%-90%)", "Small or Large\n(< 10% or > 90%)")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~zea)
# p_times
ggsave(file = png_times,
       plot = p_times,
       device = "png",
       width = 5, height = 4,
       units = "in", dpi = 300)

ggsave(file = png_times_lzw,
       plot = p_times,
       device = "tiff",
       width = 5, height = 4,
       compression = "lzw", type = "cairo",
       units = "in", dpi = 300)

# save plot as RDS also
saveRDS(object = p_times, file = rds_times)

print("summary of time estimates")
summary_times <- d %>%
  dplyr::group_by(zea, admixture_pulse) %>%
  dplyr::summarise(subspecies_min = min(time),
                 subspecies_mean = mean(time),
                 subspecies_median = median(time),
                 subspecies_max = max(time))
print(summary_times)

# write out results to a text file
left_join(d, summary_times, by = c("zea", "admixture_pulse")) %>%
write.table(., file = txt_times,
            col.names = T, row.names = F, sep = "\t", quote = F)

# so around 1300 generations on the high end.
# how many markers do I have per block?
# mean block length in cM:
mean_cM_block_length = 1/max(d$time) * 100 #cM/Morgan

allo <- read.table(paste0("local_ancestry/results/thinnedSNPs/", prefix, "/K3/whole_genome.allo.counts"),
                    stringsAsFactors = F, header = T) %>%
  dplyr::mutate(n_tot_maize = n_minor_maize + n_major_maize,
                p_maize = n_minor_maize/(n_tot_maize),
                n_tot_mex = n_minor_mex + n_major_mex,
                p_mex = n_minor_mex/(n_tot_mex),
                p_diff = abs(p_mex - p_maize))

mean_cM_between_markers = mean(allo$distM[allo$distM != 1])*100
print("Mean block length (cM):")
mean_cM_block_length
print("Mean markers per cM:")
1/mean_cM_between_markers
print("Mean markers per ancestry block:")
mean_cM_block_length/mean_cM_between_markers

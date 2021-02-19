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
prefix_all = snakemake@params[["prefix_all"]]
# prefix_all = "HILO_MAIZE55"
png_times = snakemake@output[["png_times"]]
# png_times = "local_ancestry/plots/admix_times_Ne10000_yesBoot.png"
txt_times = snakemake@output[["txt_times"]]
# txt_times = "local_ancestry/plots/admix_times_Ne10000_yesBoot.txt"
rds_times = snakemake@output[["rds_times"]]
# rds_times = "local_ancestry/plots/admix_times_Ne10000_yesBoot.txt"
source(snakemake@params[["colors"]]) # plotting colors
# source("colors.R")

# load meta and bootstrap data for maize and mexicana
zea = c("maize", "mexicana")
for (z in zea){
  load(paste0("local_ancestry/results/ancestry_hmm/",
              prefix_all, "/Ne", Ne, "_", YESNO, 
              "Boot/anc/", z, ".pop.meta.RData"))
  times = do.call(bind_cols, lapply(meta_pops$pop, function(p) 
    read.table(paste0("local_ancestry/results/ancestry_hmm/",
                      prefix_all, "/Ne", Ne, "_", YESNO, 
                      "Boot/", p, ".times"), 
               header = F, stringsAsFactors = F)))
  colnames(times) <- meta_pops$pop
  boots = times[2:nrow(times), ] # first row is estimate, remaining rows bootstraps
  time_est = data.frame(time = unlist(times[1, ]),
                    low_boot = apply(boots, 2, function(x) quantile(x, alpha/2)),
                    high_boot = apply(boots, 2, function(x) quantile(x, 1 - alpha/2)),
                    mean_boot = apply(boots, 2, mean),
                    sd_boot = apply(boots, 2, sd),
                    se_boot = sd_boot/sqrt(nrow(boots)),
                    pop = colnames(times), stringsAsFactors = F,
                    alpha = alpha) %>%
    dplyr::mutate(low_basic = 2*time - high_boot,
                  high_basic = 2*time - low_boot)
  assign(x = z, value = left_join(meta_pops, time_est, by = "pop"))
  rm(meta_pops, times, time_est, boots)
}

# combine data
d <- bind_rows(maize, mexicana) %>%
  dplyr::mutate(introgress = ifelse(zea == "maize", alpha_local_ancestry, 1 - alpha_local_ancestry))

# plot
p_times <- d %>%
  ggplot(., aes(x = reorder(LOCALITY, ELEVATION), 
                y = time,
                #alpha = introgress,
                col = zea)) +
  geom_point() + # plot time estimate
  # add errorbars for 95% percentile CI around that mean
  # based on bootstrap
  geom_errorbar(aes(ymin = low_boot,
                    ymax = high_boot),
                width = .5) +
  xlab("Sympatric population") +
  ylab("Admixture time (generations)") +
  theme_classic() +
  ylim(0, 1500) +
  guides(color = guide_legend("Subspecies"),
         alpha = guide_legend("Proportion\nminor ancestry")) +
  scale_color_manual(values = col_maize_mex_parv) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p_times
ggsave(file = png_times,
       plot = p_times,
       device = "png",
       width = 5, height = 4, 
       units = "in", dpi = 300)

# save plot as RDS also
saveRDS(object = p_times, file = rds_times)

print("summary of time estimates")
summary_times <- d %>%
  dplyr::group_by(zea) %>%
  dplyr::summarise(subspecies_min = min(time),
                 subspecies_mean = mean(time),
                 subspecies_median = median(time),
                 subspecies_max = max(time))
print(summary_times)

# write out results to a text file
left_join(d, summary_times, by = "zea") %>%
write.table(., file = txt_times, 
            col.names = T, row.names = F, sep = "\t", quote = F)

# so around 1300 generations on the high end.
# how many markers do I have per block?
# mean block length in cM:
mean_cM_block_length = 1/max(d$time) * 100 #cM/Morgan

allo <- read.table(paste0("local_ancestry/results/thinnedSNPs/", prefix_all, "/whole_genome.allo.counts"),
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

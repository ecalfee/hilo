#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
# this script plots ancestry outliers across the genome,
# from individual populations, and shared across pops

# load variables from Snakefile
#zea = snakemake@params[["zea"]]
# zea = "maize"
colors_file = snakemake@input[["colors"]]
# colors_file = "colors.R"

#sim_file = snakemake@input[["sim"]]
# sim_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".MVN.RData")
#anc_file = snakemake@input[["anc"]]
# anc_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pops.anc.RData")
#meta_file = snakemake@input[["meta_pop"]]
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

png_hist = paste0("ZAnc/plots/Ne", Ne, "_", YESNO, "Boot/mex_maize_hist_outlier_peaks.png")

# load data
source(colors_file)
# additional colors for simulated vs. observed data:
col_maize_mex_sim = c("#E69F00", col_maize_mex_parv[["maize"]],
                      col_maize_mex_parv[["mexicana"]], "plum")
names(col_maize_mex_sim) = c("maize_data", "maize_mvn_sim", 
                             "mexicana_data", "mexicana_mvn_sim")

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

# load both mexicana and maize local ancestry data
mex_maize = c("mexicana", "maize")

meta_pops_list = list(mexicana = NULL, maize = NULL)
anc_outliers_list = list(mexicana = NULL, maize = NULL)
outliers_by_pop_list = list(mexicana = NULL, maize = NULL)

for (zea in mex_maize){
  sim_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".MVN.RData")
  anc_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pops.anc.RData")
  meta_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pop.meta.RData")

  load(anc_file)
  load(sim_file)
  load(meta_file)
  
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
  
  outliers_by_pop <- bind_rows(mutate(anc_outliers, source = "data"),
                              mutate(sim_outliers, source = "mvn_sim", 
                                     inv4m = F)) %>%
    group_by(chr, pos, source, inv4m) %>%
    summarise(outlier_pops = sum(top_sd2)) %>%
    ungroup() %>%
    mutate(zea = zea)
  
  meta_pops_list[[zea]] <- meta_pops
  anc_outliers_list[[zea]] <- anc_outliers
  outliers_by_pop_list[[zea]] <- outliers_by_pop
  }

# plot local ancestry by population across each chromosome individually
for (zea in mex_maize){
  # make a plot of outliers on each chromosome
  png_chr_i_prefix = paste0("ZAnc/plots/Ne", Ne, "_", YESNO, "Boot/", zea, "_shared_outliers_chr_")
  for (i in 1:10){
    p_chr_i = anc_outliers_list[[zea]] %>%
      filter(chr == i) %>%
      mutate(top_sd2 = ifelse(top_sd2, "outlier", "non-outlier")) %>%
      # plot population ancestry at each SNP by position on chromosome (Mb)
      ggplot(., aes(pos/10^6, anc, color = top_sd2)) + 
      geom_point(size = 0.1,
                 shape = 20) +
      facet_grid(reorder(LOCALITY, desc(ELEVATION)) ~ .) +
      # add population mean introgressed ancestry
      geom_hline(data = meta_pops_list[[zea]], 
                 aes(yintercept = alpha_local_ancestry,
                     color = zea), 
                 linetype = "dashed",
                 alpha = 1
      ) +
      # add vertical lines for known inversions
      geom_vline(data = filter(inv, chr == i) %>%
                   left_join(., mutate(meta_pops_list[[zea]], chr = i), 
                             by = "chr") %>%
                   pivot_longer(cols = c("start", "end"),
                                names_to = "which_end",
                                values_to = "inv_pos"),
                 size = 0.25,
                 alpha = 1,
                 mapping = aes(xintercept = inv_pos/10^6)) +
      xlab(paste0("Chr ", i, " position (Mbp)")) +
      ylab(paste0("Introgressed ", mex_maize[mex_maize != zea], " ancestry frequency")) +
      scale_color_manual(values = c(col_maize_mex_parv[[zea]], "#00BFC4", "darkgrey"),
                         labels = c("mean ancestry", "> 2 s.d. above mean", "non-outlier"),
                         limits = c(zea, "outlier", "non-outlier")) + 
      theme_light() +
      scale_y_continuous(breaks = c(0, 1)) +
      scale_x_continuous(expand = expand_scale(mult = c(0, 0), # remove extra expanded scale 
                                               add = c(2, 2))) +
      coord_cartesian(ylim = c(0, 1)) + 
      theme(strip.text.y = element_text(angle=0),
            legend.position = "bottom",
            legend.title = element_blank(),
            panel.border = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            legend.margin = margin(c(0,1,1,0)), # top, right, bottom, left
            legend.box.margin = margin(c(-5,0,0,0))
      ) + 
      guides(color = guide_legend(override.aes = list(shape = 15, 
                                                      linetype = 0,
                                                      size = 3)))
    ggsave(file = paste0(png_chr_i_prefix, i, ".png"),
           plot = p_chr_i,
           height = 5, 
           width = 7.5, 
           units = "in",
           device = "png")
  }
}


# plot sfs of shared high-introgression outliers across pops
# using geom_bar()
p_hist <- do.call(bind_rows, outliers_by_pop_list) %>%
  group_by(outlier_pops, zea, source) %>%
  summarise(n = n()) %>%
  left_join(., 
            group_by(., zea, source) %>%
              summarise(tot = sum(n)),
            by = c("zea", "source")) %>%
  ungroup() %>%
  mutate(freq = n/tot) %>% # proportion of total snps from that zea and source
  # e.g. maize simulated data, that have that number of outlier pops
  complete(outlier_pops, zea, source,
           fill = list(freq = 0)) %>% # any categories not found have freq = 0
  filter(outlier_pops >= 1) %>% # don't plot SNPs with zero outliers
  #View()
  ggplot(., aes(x = outlier_pops, y = freq)) + #stat(width*density))) + # width*density turns it into a proportion
  geom_bar(aes(fill = paste(zea, source, sep = "_")), 
           position = "dodge2",
           stat = "identity") + 
  scale_fill_manual(values = col_maize_mex_sim,
                    name = NULL, 
                    labels = c("maize data", "maize simulations",
                               "mexicana data", "mexicana simulations")) +
  theme_light() +
  xlab("Number of populations with high introgression") +
  ylab("Proportion of SNPs") +
  coord_cartesian(xlim = c(1, 14), 
                  ylim = c(0, 0.17)) +
  facet_wrap(~zea)
ggsave(file = png_hist,
       plot = p_hist,
       height = 4, width = 6, 
       units = "in",
       device = "png") 

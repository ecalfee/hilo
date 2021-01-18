#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(tidyr)
library(widyr) # to count peaks pairwise
#library(GGally) # ggnet2
#library(network) # for graph/networks
#library(sna) # for graph/networks
# tidyverse networks/graphs and plotting:
library(tidygraph)
library(ggraph)
library(geodist) # to calculate distances between locations by lat/lon

# this script plots ancestry outliers across the genome,
# from individual populations, and shared across pops

# load variables from Snakefile
colors_file = snakemake@input[["colors"]]
# colors_file = "colors.R"
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
sim_outliers_list = list(mexicana = NULL, maize = NULL)
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
  sim_outliers_list[[zea]] <- sim_outliers
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

# ------------------------------------------------------------#
# get pairwise peak sharing
# pairwise_count(dat, letter, group, sort = TRUE, upper = F)
# (I'll need this for the simulated data too)
# now plot as network

# tidygraph/ggraph way:
nodes <- meta_pops_list[["maize"]] %>%
  arrange(ELEVATION) %>%
  mutate(id = 1:14) %>% # must be contiguous id's for nodes
  #dplyr::select(., id, pop, LOCALITY, ELEVATION) %>%
  rename(name = LOCALITY)
total_snps_data <- nrow(sites)
edges_data <- anc_outliers_list[["maize"]] %>%
  dplyr::filter(top_sd2) %>%
  dplyr::mutate(snp = paste(chr, pos, sep = "_")) %>%
  widyr::pairwise_count(., LOCALITY, snp, sort = F, upper = F) %>% # upper = F means no repeat counts for pairs
  dplyr::rename(from = item1, to = item2, n_peaks_shared_data = n) %>%
  dplyr::mutate(p_snps_shared_data = n_peaks_shared_data/total_snps_data) %>%
  arrange(from) # arrange by id

# get geographic distances between populations!

net_tidy_data <- tbl_graph(nodes = nodes, 
                      edges = edges_data, # tbl_graph() automatically substitutes node id #'s for names on each edge
                      directed = F)

# circle plot
p_net_circular <- ggraph(net_tidy_data, layout = "circle") +
  geom_edge_link(aes(width = p_snps_shared_data*100), alpha = 0.5) +
  geom_node_point(aes(color = ELEVATION), size = 4) +
  theme_graph() +
  geom_node_text(aes(label = name), repel = F) +
  scale_edge_width(range = c(0.01, 2)) +
  scale_color_viridis(direction = -1) +
  coord_cartesian(clip = "off") +
  labs(edge_width = "% SNPs shared\n introgression peaks") +
  ggtitle("Shared introgression peaks in maize")
p_net_circular
ggsave(file = "ZAnc/plots/network_maize_peak_sharing_circular.png",
       plot = p_net_circular,
       height = 4, width = 6, 
       units = "in",
       device = "png") 
# geom_node_text(nudge_x = p$data$x * .1, nudge_y = p$data$y * .1)

# by lat/lon
p_net_map <- layout_tbl_graph_manual(net_tidy_data,
                           x = LAT, y = LONG) %>%
  ggraph(.) +
  geom_edge_link(aes(width = p_snps_shared_data*100), alpha = 0.5) +
  geom_node_point(aes(color = ELEVATION), size = 4) +
  theme_graph() +
  geom_node_text(aes(label = name), repel = F) +
  scale_edge_width(range = c(0.01, 2)) +
  scale_color_viridis(direction = -1) +
  coord_cartesian(clip = "off") +
  labs(edge_width = "% SNPs shared\n introgression peaks") +
  ggtitle("Shared introgression peaks in maize")
p_net_map

# linear plot
p_net_linear <- ggraph(net_tidy_data, layout = "linear") +
  geom_edge_arc(aes(width = p_snps_shared_data*100), alpha = 0.5) +
  geom_node_point(aes(color = ELEVATION), size = 3) + # why can't I do x = ELEVATION?
  theme_graph() +
  geom_node_text(aes(label = name), 
                 angle = 90,
                 hjust = 1,
                 y = -0.5,
                 repel = F) +
  scale_edge_width(range = c(0.01, 2)) +
  scale_color_viridis(direction = -1) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(c(t = 20, r = 5.5, b = 100, l = 5.5), 
                             unit = "pt")) +
  labs(edge_width = "% SNPs shared\n introgression peaks") +
  ggtitle("Shared introgression peaks in maize")
p_net_linear
# I could color the arcs by geographic distance between pops
ggsave(file = "ZAnc/plots/network_maize_peak_sharing_linear.png",
       plot = p_net_linear,
       height = 4, width = 6, 
       units = "in",
       device = "png") 

# make linear network plot of enrichment over simulations 
# for pairwise sharing
# tidygraph/ggraph way:
total_snps_sim <- length(unique(sim_outliers_list[["maize"]]$pos))
edges_sim <- sim_outliers_list[["maize"]] %>%
  dplyr::filter(top_sd2) %>%
  dplyr::mutate(snp = paste("sim", pos, sep = "_")) %>%
  widyr::pairwise_count(., LOCALITY, snp, sort = F, upper = F) %>% # upper = F means no repeat counts for pairs
  dplyr::rename(from = item1, to = item2, n_peaks_shared_sim = n) %>%
  # proportion of snps that are shared as peaks > 2sd above mean introgression
  dplyr::mutate(p_snps_shared_sim = n_peaks_shared_sim/total_snps_sim) %>%
  arrange(from) # arrange by id
edges_both <- full_join(edges_data, edges_sim, by = c("from", "to")) %>%
  dplyr::mutate(enrichment_shared_peaks = p_snps_shared_data-p_snps_shared_sim)

net_tidy_both <- tbl_graph(nodes = nodes, 
                      edges = edges_both, # tbl_graph() automatically substitutes node id #'s for names on each edge
                      directed = F)
# linear plot
p_net_linear_enrich <- ggraph(net_tidy_both, layout = "linear") +
  geom_edge_arc(aes(width = enrichment_shared_peaks*100), alpha = 0.5) +
  geom_node_point(aes(color = ELEVATION), size = 3) + # why can't I do x = ELEVATION?
  theme_graph() +
  geom_node_text(aes(label = name), 
                 angle = 90,
                 hjust = 1,
                 y = -0.5,
                 repel = F) +
  scale_edge_width(range = c(0.01, 2)) +
  scale_color_viridis(direction = -1) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(c(t = 20, r = 5.5, b = 100, l = 5.5), 
                             unit = "pt")) +
  labs(edge_width = "Enrichment of shared peaks") +
  ggtitle("Shared introgression peaks in maize")
p_net_linear_enrich
# I could color the arcs by geographic distance between pops
ggsave(file = "ZAnc/plots/network_maize_peak_sharing_enrichment_linear.png",
       plot = p_net_linear_enrich,
       height = 4, width = 6, 
       units = "in",
       device = "png") 

# simulation data only -- linear plot
net_tidy_sim <- tbl_graph(nodes = nodes, 
                           edges = edges_sim, # tbl_graph() automatically substitutes node id #'s for names on each edge
                           directed = F)
p_net_linear_sim <- ggraph(net_tidy_sim, layout = "linear") +
  geom_edge_arc(aes(width = p_snps_shared_sim*100), alpha = 0.5) +
  geom_node_point(aes(color = ELEVATION), size = 3) + # why can't I do x = ELEVATION?
  theme_graph() +
  geom_node_text(aes(label = name), 
                 angle = 90,
                 hjust = 1,
                 y = -0.5,
                 repel = F) +
  scale_edge_width(range = c(0.01, 2)) +
  scale_color_viridis(direction = -1) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(c(t = 20, r = 5.5, b = 100, l = 5.5), 
                             unit = "pt")) +
  labs(edge_width = "% SNPs shared as peaks") +
  ggtitle("MVN simulated introgression in maize")
p_net_linear_sim
ggsave(file = "ZAnc/plots/network_maize_peak_sharing_simulations_linear.png",
       plot = p_net_linear_sim,
       height = 4, width = 6, 
       units = "in",
       device = "png") 

# get geographic distances between populations!


# --------------------------------------------------------------#
# repeat analysis defining outliers as > 50% mexicana ancestry
edges_data50 <- anc_outliers_list[["maize"]] %>%
  dplyr::filter(anc > 0.5) %>%
  dplyr::mutate(snp = paste(chr, pos, sep = "_")) %>%
  widyr::pairwise_count(., LOCALITY, snp, sort = F, upper = F) %>% # upper = F means no repeat counts for pairs
  dplyr::rename(from = item1, to = item2, n_peaks_shared_data = n) %>%
  dplyr::mutate(p_snps_shared_data = n_peaks_shared_data/total_snps_data) %>%
  arrange(from) # arrange by id

net_tidy_data50 <- tbl_graph(nodes = nodes, 
                           edges = edges_data50, # tbl_graph() automatically substitutes node id #'s for names on each edge
                           directed = F)
edges_sim50 <- sim_outliers_list[["maize"]] %>%
  dplyr::filter(anc > 0.5) %>%
  dplyr::mutate(snp = paste("sim", pos, sep = "_")) %>%
  widyr::pairwise_count(., LOCALITY, snp, sort = F, upper = F) %>% # upper = F means no repeat counts for pairs
  dplyr::rename(from = item1, to = item2, n_peaks_shared_sim = n) %>%
  # proportion of snps that are shared as peaks > 2sd above mean introgression
  dplyr::mutate(p_snps_shared_sim = n_peaks_shared_sim/total_snps_sim) %>%
  arrange(from) # arrange by id
edges_both50 <- full_join(edges_data50, edges_sim50, by = c("from", "to")) %>%
  dplyr::mutate(enrichment_shared_peaks = p_snps_shared_data/p_snps_shared_sim)

net_tidy_both50 <- tbl_graph(nodes = nodes, 
                           edges = edges_both50, # tbl_graph() automatically substitutes node id #'s for names on each edge
                           directed = F)
# linear plot of data
p_net_linear50 <- ggraph(net_tidy_data50, layout = "linear") +
  geom_edge_arc(aes(width = p_snps_shared_data*100), alpha = 0.5) +
  geom_node_point(aes(color = ELEVATION), size = 3) + # why can't I do x = ELEVATION?
  theme_graph() +
  geom_node_text(aes(label = name), 
                 angle = 90,
                 hjust = 1,
                 y = -0.5,
                 repel = F) +
  scale_edge_width(range = c(0.01, 2)) +
  scale_color_viridis(direction = -1) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(c(t = 20, r = 5.5, b = 100, l = 5.5), 
                             unit = "pt")) +
  labs(edge_width = "% SNPs shared\n introgression peaks") +
  ggtitle("Shared introgression peaks in maize")
p_net_linear50
# linear plot of enrichment
p_net_linear_enrich50 <- ggraph(net_tidy_both50, layout = "linear") +
  geom_edge_arc(aes(width = enrichment_shared_peaks), alpha = 0.5) +
  geom_node_point(aes(color = ELEVATION), size = 3) + # why can't I do x = ELEVATION?
  theme_graph() +
  geom_node_text(aes(label = name), 
                 angle = 90,
                 hjust = 1,
                 y = -0.5,
                 repel = F) +
  scale_edge_width(range = c(0.01, 2)) +
  scale_color_viridis(direction = -1) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(c(t = 20, r = 5.5, b = 100, l = 5.5), 
                             unit = "pt")) +
  labs(edge_width = "Enrichment of shared peaks") +
  ggtitle("Shared introgression peaks in maize")
p_net_linear_enrich50

# simulation data only -- linear plot
net_tidy_sim50 <- tbl_graph(nodes = nodes, 
                          edges = edges_sim50, # tbl_graph() automatically substitutes node id #'s for names on each edge
                          directed = F)
p_net_linear_sim50 <- ggraph(net_tidy_sim50, layout = "linear") +
  geom_edge_arc(aes(width = p_snps_shared_sim*100), alpha = 0.5) +
  geom_node_point(aes(color = ELEVATION), size = 3) + # why can't I do x = ELEVATION?
  theme_graph() +
  geom_node_text(aes(label = name), 
                 angle = 90,
                 hjust = 1,
                 y = -0.5,
                 repel = F) +
  scale_edge_width(range = c(0.01, 2)) +
  scale_color_viridis(direction = -1) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(c(t = 20, r = 5.5, b = 100, l = 5.5), 
                             unit = "pt")) +
  labs(edge_width = "% SNPs shared as peaks") +
  ggtitle("MVN simulated introgression in maize")
p_net_linear_sim50

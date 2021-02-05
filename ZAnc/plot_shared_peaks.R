#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(tidyr)
# to count peaks pairwise:
library(widyr)
# tidyverse networks/graphs and plotting:
library(tidygraph)
library(ggraph)
# to calculate distances between locations by lat/lon:
library(geodist) 
# for multi-panel plots:
library(grid)
library(gridExtra)
library(cowplot)
library(forcats)

# this script plots ancestry outliers across the genome,
# from individual populations, and shared across pops

# load variables from Snakefile
colors_file = snakemake@input[["colors"]]
# colors_file = "colors.R"
sites_file = snakemake@input[["sites"]]
# sites_file = "local_ancestry/results/thinnedSNPs/HILO_MAIZE55/whole_genome.var.sites"
inv_file = snakemake@input[["inv"]]
# inv_file = "data/refMaize/inversions/knownInv_v4_coord.txt"
genome_file = snakemake@input[["genome"]]
# genome_file = "data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"
Ne = as.integer(snakemake@params[["Ne"]])
# Ne = 10000
YESNO = snakemake@params[["YESNO"]]
# YESNO = "yes"
prefix_all = snakemake@params[["prefix_all"]]
# prefix_all = "HILO_MAIZE55"

png_hist = paste0("ZAnc/plots/Ne", Ne, "_", YESNO, "Boot/mex_maize_hist_outlier_peaks.png")
png_net_multi = paste0("ZAnc/plots/Ne", Ne, "_", YESNO, "Boot/network_peak_sharing.png")
png_net_multi_no_inv4m = paste0("ZAnc/plots/Ne", Ne, "_", YESNO, "Boot/network_peak_sharing_no_inv4m.png")
png_net_multi_no_inv4m_data = paste0("ZAnc/plots/Ne", Ne, "_", YESNO, "Boot/network_peak_sharing_no_inv4m_data_only.png")
png_net_multi_sims = paste0("ZAnc/plots/Ne", Ne, "_", YESNO, "Boot/network_peak_sharing_sims_only.png")
png_combmatrix_maize = paste0("ZAnc/plots/Ne", Ne, "_", YESNO, "Boot/combmatrix_peak_sharing_maize.png")
png_combmatrix_mexicana = paste0("ZAnc/plots/Ne", Ne, "_", YESNO, "Boot/combmatrix_peak_sharing_mexicana.png")
png_peaks_on_genome = paste0("ZAnc/plots/Ne", Ne, "_", YESNO, "Boot/shared_peaks_on_genome.png")

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

# load chromosome lengths
genome <- read.table(genome_file, header = F, stringsAsFactors = F,
                     sep = "\t") %>%
  data.table::setnames(c("chr", "length")) %>%
  dplyr::mutate(chr_end = cumsum(length),
                chr_start = c(0, chr_end[1:(nrow(.)-1)]))

# and site/position information for SNPs with ancestry calls
sites <- read.table(sites_file, header = F, stringsAsFactors = F,
                    sep = "\t") %>%
  data.table::setnames(c("chr", "pos", "major", "minor")) %>%
  dplyr::mutate(inv4m = (chr == inv$chr[inv$ID=="inv4m"] & 
                           pos >= inv$start[inv$ID=="inv4m"] & 
                           pos <= inv$end[inv$ID == "inv4m"])) %>% # mark sites within inversion
  left_join(., genome, by = "chr") %>%
  dplyr::mutate(pos_cum = chr_start + pos) # get cumulative chromosomal position

# mean of each chromosome (cumulative position)
axis_spacing = sites %>%
  group_by(chr) %>% 
  summarize(center=(max(pos_cum) + min(pos_cum)) / 2,
            start = min(pos_cum),
            end = max(pos_cum))

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
# p_hist
ggsave(file = png_hist,
       plot = p_hist,
       height = 4, width = 6, 
       units = "in",
       device = "png") 

# ------------------------------------------------------------#
# get pairwise peak sharing and plot as a network
peak_network = list(mexicana = NULL, maize = NULL)

for (zea in mex_maize){

  nodes <- meta_pops_list[[zea]] %>%
    arrange(ELEVATION) %>%
    mutate(id = 1:14) %>% # must be contiguous id's for nodes
    rename(name = LOCALITY)
  
  # count shared outlier peaks in data
  total_snps_data <- nrow(sites)
  edges_data <- anc_outliers_list[[zea]] %>%
    dplyr::filter(top_sd2) %>%
    dplyr::mutate(snp = paste(chr, pos, sep = "_")) %>%
    widyr::pairwise_count(., LOCALITY, snp, sort = F, upper = F) %>% # upper = F means no repeat counts for pairs
    dplyr::rename(from = item1, to = item2, n_peaks_shared_data = n) %>%
    dplyr::mutate(p_snps_shared_data = n_peaks_shared_data/total_snps_data) %>%
    arrange(from) # arrange by id
  
  # count shared outlier peaks in simulations
  total_snps_sim <- length(unique(sim_outliers_list[[zea]]$pos))
  edges_sim <- sim_outliers_list[[zea]] %>%
    dplyr::filter(top_sd2) %>%
    dplyr::mutate(snp = paste("sim", pos, sep = "_")) %>%
    widyr::pairwise_count(., LOCALITY, snp, sort = F, upper = F) %>% # upper = F means no repeat counts for pairs
    dplyr::rename(from = item1, to = item2, n_peaks_shared_sim = n) %>%
    # proportion of snps that are shared as peaks > 2sd above mean introgression
    dplyr::mutate(p_snps_shared_sim = n_peaks_shared_sim/total_snps_sim) %>%
    arrange(from) # arrange by id
  
  # edges data compared to simulations
  edges_both <- full_join(edges_data, edges_sim, by = c("from", "to")) %>%
    # how much more sharing do we see in real data compared to simulations?
    dplyr::mutate(surplus_shared_peaks = p_snps_shared_data - p_snps_shared_sim)
  
  # create reverse network from <-> to
  edges_both_reverse <- edges_both %>%
    rename(from_old = from, to_old = to) %>%
    mutate(from = to_old, to = from_old) %>%
    dplyr::select(from, to, surplus_shared_peaks, p_snps_shared_data, p_snps_shared_sim)
    
  if (zea == "maize"){
    net_tidy <- tbl_graph(nodes = nodes, # tbl_graph() automatically substitutes node id #'s for names on each edge
                          edges = edges_both, 
                          directed = T)
  }else{
    net_tidy <- tbl_graph(nodes = nodes, 
                          edges = edges_both_reverse, # reverse direction of edges if mexicana (to plot upsidedown)
                          directed = T)
  }
  
  peak_network[[zea]] <- net_tidy
}

# linear plot
p_net_maize <- ggraph(peak_network[["maize"]], layout = "linear") +
    geom_edge_arc(aes(width = surplus_shared_peaks*100,
                      alpha = surplus_shared_peaks*100)) +
    geom_node_point(aes(color = ELEVATION), size = 3) + # why can't I do x = ELEVATION?
    theme_graph(base_family = 'Helvetica') +
    geom_node_text(aes(label = name), 
                   angle = 90,
                   hjust = 1,
                   y = -0.5,
                   repel = F) +
    scale_edge_alpha(range = c(0, 1), limits = c(0, 2.5)) +
    scale_edge_width(range = c(0, 2.5), limits = c(0, 2.5)) +
    scale_color_viridis(direction = -1) +
    coord_cartesian(clip = "off") +
  labs(edge_width = "% SNPs in\n shared peaks\n beyond expectation",
       edge_alpha = "% SNPs in\n shared peaks\n beyond expectation",
       color = "Elevation (m)") +
  ggtitle("Shared introgression peaks in maize")
# p_net_maize

p_net_mexicana <- ggraph(peak_network[["mexicana"]], layout = "linear") +
  geom_edge_arc(aes(width = surplus_shared_peaks*100, 
                    alpha = surplus_shared_peaks*100)) +
  geom_node_point(aes(color = ELEVATION), size = 3) + # why can't I do x = ELEVATION?
  theme_graph(base_family = 'Helvetica') +
  scale_edge_alpha(range = c(0, 1), limits = c(0, 2.5)) +
  scale_edge_width(range = c(0, 2.5), limits = c(0, 2.5)) +
  scale_color_viridis(direction = -1) +
  coord_cartesian(clip = "off") +
  labs(edge_width = "% SNPs in\n shared peaks\n beyond expectation",
       edge_alpha = "% SNPs in\n shared peaks\n beyond expectation",
       color = "Elevation (m)") +
  ggtitle("Shared introgression peaks in mexicana")
# p_net_mexicana

# combine maize and mexicana networks into 1 plot:
p_net_multi <- grid.arrange(grobs = list(ggplotGrob(p_net_maize +
                                                      #labs(subtitle = "maize") +
                                                      theme(plot.margin = margin(c(t = 0, r = 5, b = 95, l = 2.5), unit = "pt"),
                                                            plot.title = element_blank(),
                                                            legend.position = "none"
                                                            )),
                                     ggplotGrob(p_net_mexicana +
                                                  #labs(subtitle = "mexicana") +
                                                  theme(plot.margin = margin(c(t = 5, r = 5, b = 5, l = 2.5), unit = "pt"),
                                                        plot.title = element_blank(),
                                                        legend.position = "none"
                                                        )),
                                     cowplot::get_legend(p_net_maize),
                                     textGrob(label = "maize", 
                                              x = unit(1, "lines"), 
                                              y = unit(10, "lines"),
                                              rot = 90),
                                     textGrob(label = "mexicana", 
                                              x = unit(1, "lines"), 
                                              y = unit(9, "lines"),
                                              rot = 90)),
                            layout_matrix = rbind(c(4, 1, NA, 3),
                                                  c(5, 2, NA, 3)),
                            heights = c(1, 
                                        0.7),
                            widths = c(0.1, 5, 0.2, 1.7))

# p_net_multi

ggsave(file = png_net_multi,
         plot = p_net_multi,
         height = 6.5, width = 6.2, 
         units = "in",
         device = "png") 



# ------------------ network excluding known inv4m ---------------------
# get pairwise peak sharing and plot as a network
peak_network_no_inv4m = list(mexicana = NULL, maize = NULL)

for (zea in mex_maize){
  
  nodes <- meta_pops_list[[zea]] %>%
    arrange(ELEVATION) %>%
    mutate(id = 1:14) %>% # must be contiguous id's for nodes
    rename(name = LOCALITY)
  
  # count shared outlier peaks in data
  total_snps_data <- nrow(filter(sites, !inv4m))
  edges_data <- anc_outliers_list[[zea]] %>%
    dplyr::filter(top_sd2) %>%
    dplyr::filter(!inv4m) %>%
    dplyr::mutate(snp = paste(chr, pos, sep = "_")) %>%
    widyr::pairwise_count(., LOCALITY, snp, sort = F, upper = F) %>% # upper = F means no repeat counts for pairs
    dplyr::rename(from = item1, to = item2, n_peaks_shared_data = n) %>%
    dplyr::mutate(p_snps_shared_data = n_peaks_shared_data/total_snps_data) %>%
    arrange(from) # arrange by id
  
  # count shared outlier peaks in simulations
  total_snps_sim <- length(unique(sim_outliers_list[[zea]]$pos))
  edges_sim <- sim_outliers_list[[zea]] %>%
    dplyr::filter(top_sd2) %>%
    dplyr::mutate(snp = paste("sim", pos, sep = "_")) %>%
    widyr::pairwise_count(., LOCALITY, snp, sort = F, upper = F) %>% # upper = F means no repeat counts for pairs
    dplyr::rename(from = item1, to = item2, n_peaks_shared_sim = n) %>%
    # proportion of snps that are shared as peaks > 2sd above mean introgression
    dplyr::mutate(p_snps_shared_sim = n_peaks_shared_sim/total_snps_sim) %>%
    arrange(from) # arrange by id
  
  # edges data compared to simulations
  edges_both <- full_join(edges_data, edges_sim, by = c("from", "to")) %>%
    # how much more sharing do we see in real data compared to simulations?
    dplyr::mutate(surplus_shared_peaks = p_snps_shared_data - p_snps_shared_sim)
  
  print(zea)
  print("raw data (%):")
  print(summary(edges_both$p_snps_shared_data*100))
  print("surplus sharing (%):")
  print(summary(edges_both$surplus_shared_peaks*100))
  # create reverse network from <-> to
  edges_both_reverse <- edges_both %>%
    rename(from_old = from, to_old = to) %>%
    mutate(from = to_old, to = from_old) %>%
    dplyr::select(from, to, surplus_shared_peaks, p_snps_shared_data, p_snps_shared_sim)
  
  if (zea == "maize"){
    net_tidy <- tbl_graph(nodes = nodes, # tbl_graph() automatically substitutes node id #'s for names on each edge
                          edges = edges_both, 
                          directed = T)
  }else{
    net_tidy <- tbl_graph(nodes = nodes, 
                          edges = edges_both_reverse, # reverse direction of edges if mexicana (to plot upsidedown)
                          directed = T)
  }
  
  peak_network_no_inv4m[[zea]] <- net_tidy
}

# linear plot
p_net_maize_no_inv4m <- ggraph(peak_network_no_inv4m[["maize"]], layout = "linear") +
  geom_edge_arc(aes(width = surplus_shared_peaks*100,
                    alpha = surplus_shared_peaks*100)) +
  geom_node_point(aes(color = ELEVATION), size = 3) + # why can't I do x = ELEVATION?
  theme_graph(base_family = 'Helvetica') +
  geom_node_text(aes(label = name), 
                 angle = 90,
                 hjust = 1,
                 y = -0.5,
                 repel = F) +
  scale_edge_alpha(range = c(0, 1.5), limits = c(0, 2.75)) +
  scale_edge_width(range = c(0, 4), limits = c(0, 2.75)) +
  scale_color_viridis(direction = -1) +
  coord_cartesian(clip = "off") +
  labs(edge_width = "% SNPs in\n shared peaks\n beyond expectation",
       edge_alpha = "% SNPs in\n shared peaks\n beyond expectation",
       color = "Elevation (m)") +
  ggtitle("(Excl. inv4m) Shared introgression peaks in maize")
# p_net_maize_no_inv4m

p_net_mexicana_no_inv4m <- ggraph(peak_network_no_inv4m[["mexicana"]], layout = "linear") +
  geom_edge_arc(aes(width = surplus_shared_peaks*100, 
                    alpha = surplus_shared_peaks*100)) +
  geom_node_point(aes(color = ELEVATION), size = 3) + # why can't I do x = ELEVATION?
  theme_graph(base_family = 'Helvetica') +
  scale_edge_alpha(range = c(0, 1.5), limits = c(0, 2.75)) +
  scale_edge_width(range = c(0, 4), limits = c(0, 2.75)) +
  scale_color_viridis(direction = -1) +
  coord_cartesian(clip = "off") +
  labs(edge_width = "% SNPs in\n shared peaks\n beyond expectation",
       edge_alpha = "% SNPs in\n shared peaks\n beyond expectation",
       color = "Elevation (m)") +
  ggtitle("(Excl. inv4m) Shared introgression peaks in mexicana")
# p_net_mexicana_no_inv4m

# combine maize and mexicana networks into 1 plot:
p_net_multi_no_inv4m <- grid.arrange(grobs = list(ggplotGrob(p_net_maize_no_inv4m +
                                                      #labs(subtitle = "maize") +
                                                      theme(plot.margin = margin(c(t = 0, r = 5, b = 95, l = 2.5), unit = "pt"),
                                                            plot.title = element_blank(),
                                                            legend.position = "none"
                                                      )),
                                         ggplotGrob(p_net_mexicana_no_inv4m +
                                                      #labs(subtitle = "mexicana") +
                                                      theme(plot.margin = margin(c(t = 5, r = 5, b = 5, l = 2.5), unit = "pt"),
                                                            plot.title = element_blank(),
                                                            legend.position = "none"
                                                      )),
                                         cowplot::get_legend(p_net_maize_no_inv4m),
                                         textGrob(label = "maize", 
                                                  x = unit(1, "lines"), 
                                                  y = unit(10, "lines"),
                                                  rot = 90),
                                         textGrob(label = "mexicana", 
                                                  x = unit(1, "lines"), 
                                                  y = unit(9, "lines"),
                                                  rot = 90)),
                            layout_matrix = rbind(c(4, 1, NA, 3),
                                                  c(5, 2, NA, 3)),
                            heights = c(1, 
                                        0.7),
                            widths = c(0.1, 5, 0.2, 1.7))

# p_net_multi_no_inv4m

ggsave(file = png_net_multi_no_inv4m,
       plot = p_net_multi_no_inv4m,
       height = 6.5, width = 6.2, 
       units = "in",
       device = "png") 

# ------plot original data network excluding inv4m (not surplus, ignore simulations) ----
p_net_maize_no_inv4m_data <- ggraph(peak_network_no_inv4m[["maize"]], layout = "linear") +
  geom_edge_arc(aes(width = p_snps_shared_data*100,
                    alpha = p_snps_shared_data*100)) +
  geom_node_point(aes(color = ELEVATION), size = 3) + # why can't I do x = ELEVATION?
  theme_graph(base_family = 'Helvetica') +
  geom_node_text(aes(label = name), 
                 angle = 90,
                 hjust = 1,
                 y = -0.5,
                 repel = F) +
  scale_edge_alpha(range = c(0, 1), limits = c(0, 2.75)) +
  scale_edge_width(range = c(0, 2.5), limits = c(0, 2.75)) +
  scale_color_viridis(direction = -1) +
  coord_cartesian(clip = "off") +
  labs(edge_width = "% SNPs in\n shared peaks\n (raw data)",
       edge_alpha = "% SNPs in\n shared peaks\n (raw data)",
       color = "Elevation (m)") +
  ggtitle("(Excl. inv4m) Shared introgression peaks in maize (raw data)")
# p_net_maize_no_inv4m_data

p_net_mexicana_no_inv4m_data <- ggraph(peak_network_no_inv4m[["mexicana"]], layout = "linear") +
  geom_edge_arc(aes(width = p_snps_shared_data*100, 
                    alpha = p_snps_shared_data*100)) +
  geom_node_point(aes(color = ELEVATION), size = 3) + # why can't I do x = ELEVATION?
  theme_graph(base_family = 'Helvetica') +
  scale_edge_alpha(range = c(0, 1), limits = c(0, 2.75)) +
  scale_edge_width(range = c(0, 2.5), limits = c(0, 2.75)) +
  scale_color_viridis(direction = -1) +
  coord_cartesian(clip = "off") +
  labs(edge_width = "% SNPs in\n shared peaks\n (raw data)",
       edge_alpha = "% SNPs in\n shared peaks\n (raw data)",
       color = "Elevation (m)") +
  ggtitle("(Excl. inv4m) Shared introgression peaks in mexicana (raw data)")
# p_net_mexicana_no_inv4m_data
# combine maize and mexicana networks into 1 plot:
p_net_multi_no_inv4m_data <- grid.arrange(grobs = list(ggplotGrob(p_net_maize_no_inv4m_data +
                                                               #labs(subtitle = "maize") +
                                                               theme(plot.margin = margin(c(t = 0, r = 5, b = 95, l = 2.5), unit = "pt"),
                                                                     plot.title = element_blank(),
                                                                     legend.position = "none"
                                                               )),
                                                  ggplotGrob(p_net_mexicana_no_inv4m_data +
                                                               #labs(subtitle = "mexicana") +
                                                               theme(plot.margin = margin(c(t = 5, r = 5, b = 5, l = 2.5), unit = "pt"),
                                                                     plot.title = element_blank(),
                                                                     legend.position = "none"
                                                               )),
                                                  cowplot::get_legend(p_net_maize_no_inv4m_data),
                                                  textGrob(label = "maize", 
                                                           x = unit(1, "lines"), 
                                                           y = unit(10, "lines"),
                                                           rot = 90),
                                                  textGrob(label = "mexicana", 
                                                           x = unit(1, "lines"), 
                                                           y = unit(9, "lines"),
                                                           rot = 90)),
                                     layout_matrix = rbind(c(4, 1, NA, 3),
                                                           c(5, 2, NA, 3)),
                                     heights = c(1, 
                                                 0.7),
                                     widths = c(0.1, 5, 0.2, 1.7))

# p_net_multi_no_inv4m_data

ggsave(file = png_net_multi_no_inv4m_data,
       plot = p_net_multi_no_inv4m_data,
       height = 6.5, width = 6.2, 
       units = "in",
       device = "png") 

# ------plot simulated network data only ------------------------
p_net_maize_sims <- ggraph(peak_network[["maize"]], layout = "linear") +
  geom_edge_arc(aes(width = p_snps_shared_sim*100,
                    alpha = p_snps_shared_sim*100)) +
  geom_node_point(aes(color = ELEVATION), size = 3) + # why can't I do x = ELEVATION?
  theme_graph(base_family = 'Helvetica') +
  geom_node_text(aes(label = name), 
                 angle = 90,
                 hjust = 1,
                 y = -0.5,
                 repel = F) +
  scale_edge_alpha(range = c(0, 1), limits = c(0, 2.5)) +
  scale_edge_width(range = c(0, 2.5), limits = c(0, 2.5)) +
  scale_color_viridis(direction = -1) +
  coord_cartesian(clip = "off") +
  labs(edge_width = "% SNPs in\n shared peaks\n (simulations)",
       edge_alpha = "% SNPs in\n shared peaks\n (simulations)",
       color = "Elevation (m)") +
  ggtitle("(Simulations) Shared introgression peaks in maize")
# p_net_maize_sims

p_net_mexicana_sims <- ggraph(peak_network[["mexicana"]], layout = "linear") +
  geom_edge_arc(aes(width = p_snps_shared_sim*100, 
                    alpha = p_snps_shared_sim*100)) +
  geom_node_point(aes(color = ELEVATION), size = 3) + # why can't I do x = ELEVATION?
  theme_graph(base_family = 'Helvetica') +
  scale_edge_alpha(range = c(0, 1), limits = c(0, 2.5)) +
  scale_edge_width(range = c(0, 2.5), limits = c(0, 2.5)) +
  scale_color_viridis(direction = -1) +
  coord_cartesian(clip = "off") +
  labs(edge_width = "% SNPs in\n shared peaks\n simulations",
       edge_alpha = "% SNPs in\n shared peaks\n simulations",
       color = "Elevation (m)") +
  ggtitle("(Simulations) Shared introgression peaks in mexicana")
# p_net_mexicana_sims
# combine maize and mexicana networks into 1 plot:
p_net_multi_sims <- grid.arrange(grobs = list(ggplotGrob(p_net_maize_sims +
                                                                    #labs(subtitle = "maize") +
                                                                    theme(plot.margin = margin(c(t = 0, r = 5, b = 95, l = 2.5), unit = "pt"),
                                                                          plot.title = element_blank(),
                                                                          legend.position = "none"
                                                                    )),
                                                       ggplotGrob(p_net_mexicana_sims +
                                                                    #labs(subtitle = "mexicana") +
                                                                    theme(plot.margin = margin(c(t = 5, r = 5, b = 5, l = 2.5), unit = "pt"),
                                                                          plot.title = element_blank(),
                                                                          legend.position = "none"
                                                                    )),
                                                       cowplot::get_legend(p_net_maize_sims),
                                                       textGrob(label = "maize", 
                                                                x = unit(1, "lines"), 
                                                                y = unit(10, "lines"),
                                                                rot = 90),
                                                       textGrob(label = "mexicana", 
                                                                x = unit(1, "lines"), 
                                                                y = unit(9, "lines"),
                                                                rot = 90)),
                                          layout_matrix = rbind(c(4, 1, NA, 3),
                                                                c(5, 2, NA, 3)),
                                          heights = c(1, 
                                                      0.7),
                                          widths = c(0.1, 5, 0.2, 1.7))

# p_net_multi_sims

ggsave(file = png_net_multi_sims,
       plot = p_net_multi_sims,
       height = 6.5, width = 6.2, 
       units = "in",
       device = "png") 


# ----------which pops commonly share peaks? plot ggupset or other combination matrix ---------
maize_pops_shared <- anc_outliers_list[["maize"]] %>% 
  arrange(ELEVATION) %>%
  mutate(LOCALITY_ELEVATION = paste0(LOCALITY,
                                     " ",
                                     ELEVATION,
                                     "m"),
         LOCALITY = forcats::fct_reorder(as.factor(LOCALITY), ELEVATION)) %>%
  filter(top_sd2) %>% 
  group_by(chr, pos) %>% 
  summarize(populations = paste0(sort(LOCALITY), collapse = "-"),
            n = length(unique(LOCALITY)))

p_combmatrix_maize <- maize_pops_shared %>%
  group_by(populations, n) %>%
  summarise(freq = n()/nrow(sites)*100) %>%
  ungroup() %>%
  mutate(populations = forcats::fct_reorder(as.factor(populations), -freq)) %>%
  arrange(populations) %>%
  head(n = 100) %>%
  ggplot(data = ., mapping = aes(x = populations, y = freq, fill = n)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  axis_combmatrix(sep = "-", 
                  levels = meta_pops_list[["maize"]] %>%
                    arrange(desc(ELEVATION)) %>%
                    .$LOCALITY) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme_combmatrix(combmatrix.panel.point.color.fill = "#00BFC4",
                   combmatrix.panel.point.color.empty = "darkgrey",
                   combmatrix.panel.line.size = 0,
                   combmatrix.panel.point.size = 1.5,
                   combmatrix.label.make_space = FALSE) +
  scale_fill_viridis_c(option = "magma",
                       limits = c(1,14)) + 
  theme(plot.margin = margin(t = 5, b = 5, l = 50, r = 5, unit = "pt"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  xlab("populations sharing high introgression peaks") +
  ylab("% SNPs within peaks") +
  labs(fill = "number of\npopulations")
# p_combmatrix_maize
ggsave(file = png_combmatrix_maize,
       plot = p_combmatrix_maize,
       height = 5, width = 7.5, 
       units = "in",
       device = "png") 

# now mexicana
mexicana_pops_shared <- anc_outliers_list[["mexicana"]] %>% 
  arrange(ELEVATION) %>%
  mutate(LOCALITY_ELEVATION = paste0(LOCALITY,
                                     " ",
                                     ELEVATION,
                                     "m"),
         LOCALITY = forcats::fct_reorder(as.factor(LOCALITY), ELEVATION)) %>%
  filter(top_sd2) %>% 
  group_by(chr, pos) %>% 
  summarize(populations = paste0(sort(LOCALITY), collapse = "-"),
            n = length(unique(LOCALITY)))

p_combmatrix_mexicana <- mexicana_pops_shared %>%
  group_by(populations, n) %>%
  summarise(freq = n()/nrow(sites)*100) %>%
  ungroup() %>%
  mutate(populations = forcats::fct_reorder(as.factor(populations), -freq)) %>%
  arrange(populations) %>%
  head(n = 100) %>%
  ggplot(data = ., mapping = aes(x = populations, y = freq, fill = n)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  axis_combmatrix(sep = "-", 
                  levels = meta_pops_list[["mexicana"]] %>%
                    arrange(desc(ELEVATION)) %>%
                    .$LOCALITY) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme_combmatrix(combmatrix.panel.point.color.fill = "#00BFC4",
                   combmatrix.panel.point.color.empty = "darkgrey",
                   combmatrix.panel.line.size = 0,
                   combmatrix.panel.point.size = 1.5,
                   combmatrix.label.make_space = FALSE) +
  scale_fill_viridis_c(option = "magma",
                       limits = c(1,14)) +
  theme(plot.margin = margin(t = 5, b = 5, l = 50, r = 5, unit = "pt"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  xlab("populations sharing high introgression peaks") +
  ylab("% SNPs within peaks") +
  labs(fill = "number of\npopulations")
# p_combmatrix_mexicana
ggsave(file = png_combmatrix_mexicana,
       plot = p_combmatrix_mexicana,
       height = 5, width = 7.5, 
       units = "in",
       device = "png") 

# ----------plot where shared peaks are across the genome-------
# defined by 2 sd > mean but also maybe a cutoff like 50% introgressed ancestry
# or maybe top 1% of the genome for that population (?)


peaks_on_genome <- left_join(sites, 
                             bind_rows(outliers_by_pop_list[["maize"]],
                                       outliers_by_pop_list[["mexicana"]]),
          by = c("chr", "pos", "inv4m"))

peaks_on_genome %>%
  ggplot(., aes(x = pos_cum, y = outlier_pops, color = outlier_pops)) +
  geom_point(size = 0.1) +
  scale_x_continuous(label = axis_spacing$chr, 
                   breaks = axis_spacing$center,
                   expand = expansion(mult = c(0, 0),
                                      add = c(0, 0))) +
  theme_classic() +
  facet_grid(zea~.) +
  scale_color_viridis_c(limits = c(1, 14), option = "magma") +
  labs(color = "number of\npopulations") +
  xlab("position on genome") +
  ylab("number of populations")

top_shared_pops <- maize_pops_shared %>%
  group_by(populations, n) %>%
  summarise(freq = n()/nrow(sites)*100) %>%
  ungroup() %>%
  mutate(populations = forcats::fct_reorder(as.factor(populations), -freq)) %>%
  arrange(populations) %>%
  head(n = 100) %>%
  dplyr::select(populations, freq) 

top_shared_pops %>%
  left_join(., maize_pops_shared, by = "populations") %>%
  filter(n >= 5) %>%
  mutate(populations = forcats::fct_reorder(as.factor(populations), freq)) %>%
  inner_join(sites, ., by = c("chr", "pos")) %>%
  ggplot(., aes(x = pos_cum, y = populations, color = n)) +
  geom_point(size = 0.1) +
  scale_x_continuous(label = axis_spacing$chr, 
                     breaks = axis_spacing$center,
                     expand = expansion(mult = c(0, 0),
                                        add = c(0, 0))) +
  scale_color_viridis_c(limits = c(1, 14), option = "magma") +
  labs(color = "number of\npopulations") +
  xlab("position on genome") +
  ylab("populations sharing peak") +
  theme_light()
# this shows it's not just the inversions :) although they definitely contribute to the excess

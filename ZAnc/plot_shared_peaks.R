#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(tidyr)
# to count peaks pairwise:
library(widyr)
# tidyverse networks/graphs and plotting:
library(tidygraph)
library(ggraph)
# for multi-panel plots:
library(grid)
library(gridExtra)
library(cowplot)
library(ggupset)

# this script plots ancestry outliers across the genome,
# from individual populations, and shared across pops

# load variables from Snakefile
colors_file = snakemake@input[["colors"]]
# colors_file = "colors.R"
sites_file = snakemake@input[["sites"]]
# sites_file = "local_ancestry/results/thinnedSNPs/HILO_MAIZE55_PARV50/K3/whole_genome.var.sites"
inv_file = snakemake@input[["inv"]]
# inv_file = "data/refMaize/inversions/knownInv_v4_coord.txt"
genome_file = snakemake@input[["genome"]]
# genome_file = "data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"
Ne = as.integer(snakemake@params[["Ne"]])
# Ne = 10000
YESNO = snakemake@params[["YESNO"]]
# YESNO = "yes"
prefix = snakemake@params[["PREFIX"]]
# prefix = "HILO_MAIZE55_PARV50"
K = snakemake@params[["K"]]
# K = 3

png_net_multi_data = paste0("ZAnc/plots/", prefix, "_K", K, "_Ne", Ne, "_", YESNO, "Boot_network_peak_sharing_data_only.png")
png_net_multi_data_lzw = paste0("../hilo_manuscript/figures_main/", prefix, "_K", K, "_Ne", Ne, "_", YESNO, "Boot_network_peak_sharing_data_only.tif")

png_combmatrix_maize = paste0("ZAnc/plots/", prefix, "_K", K, "_Ne", Ne, "_", YESNO, "Boot_combmatrix_peak_sharing_maize.png")
png_combmatrix_maize_lzw = paste0("../hilo_manuscript/figures_supp/", prefix, "_K", K, "_Ne", Ne, "_", YESNO, "Boot_combmatrix_peak_sharing_maize.tif")

png_combmatrix_mexicana = paste0("ZAnc/plots/", prefix, "_K", K, "_Ne", Ne, "_", YESNO, "Boot_combmatrix_peak_sharing_mexicana.png")
png_combmatrix_mexicana_lzw = paste0("../hilo_manuscript/figures_supp/", prefix, "_K", K, "_Ne", Ne, "_", YESNO, "Boot_combmatrix_peak_sharing_mexicana.tif")

png_peaks_on_genome = paste0("ZAnc/plots/", prefix, "_K", K, "_Ne", Ne, "_", YESNO, "Boot_shared_peaks_on_genome.png")

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
  sim_file = paste0("ZAnc/results/", prefix, "/K", K, "/Ne", Ne, "_", YESNO, "Boot/", zea, ".MVN.RData")
  anc_file = paste0("local_ancestry/results/ancestry_hmm/", prefix, "/K", K, "/Ne", Ne, "_", YESNO, "Boot/anc/", zea, ".pops.anc.RData")
  meta_file = paste0("local_ancestry/results/ancestry_hmm/", prefix, "/K", K, "/Ne", Ne, "_", YESNO, "Boot/anc/", zea, ".pop.meta.RData")

  load(anc_file)
  load(sim_file)
  load(meta_file)
  
  # take inverse if mexicana to get 'introgressed' minor ancestry
  if (zea == "mexicana"){ # sympatric mexicana look at introgressed maize
    anc = anc[["maize"]]
    MVN_sim = MVN_sim[["maize"]]
    meta_pops$alpha_local_ancestry = meta_pops$alpha_local_ancestry_maize
  } else{ # sympatric maize look at introgressed mexicana
    anc = anc[["mexicana"]]
    MVN_sim = MVN_sim[["mexicana"]]
    meta_pops$alpha_local_ancestry = meta_pops$alpha_local_ancestry_mexicana
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
  png_chr_i_prefix = paste0("ZAnc/plots/", prefix, "_K", K, "_Ne", Ne, "_", YESNO, "Boot_", zea, "_shared_outliers_chr_")
  png_chr_i_prefix_main_lzw = paste0("../hilo_manuscript/figures_main/", prefix, "_K", K, "_Ne", Ne, "_", YESNO, "Boot_", zea, "_shared_outliers_chr_")
  png_chr_i_prefix_supp_lzw = paste0("../hilo_manuscript/figures_supp/", prefix, "_K", K, "_Ne", Ne, "_", YESNO, "Boot_", zea, "_shared_outliers_chr_")
  
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
            legend.key.width = unit(1.2,"cm"),
            legend.title = element_blank(),
            panel.border = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            legend.margin = margin(c(0,1,1,0)), # top, right, bottom, left
            legend.box.margin = margin(c(-5,0,0,0))
      ) + 
      guides(color = guide_legend(override.aes = list(shape = c(NA, 15, 15), 
                                                      linetype = c(2, 0, 0),
                                                      size = c(1, 3, 3))))
    if (i == 4){
      p_chr_i = p_chr_i +
        # add vertical lines for known inversion at inv4m
        geom_vline(data = filter(inv, chr == i) %>%
                     left_join(., mutate(meta_pops_list[[zea]], chr = i), 
                               by = "chr") %>%
                     pivot_longer(cols = c("start", "end"),
                                  names_to = "which_end",
                                  values_to = "inv_pos"),
                   size = 0.25,
                   alpha = 1,
                   mapping = aes(xintercept = inv_pos/10^6))
    }

    p_chr_i_arrow = grid.arrange(grobs = list(ggplotGrob(p_chr_i),
                                       textGrob(label = expression("Low to high elevation " %->% ""), 
                                                rot = 90,
                                                just = "left",
                                                gp = gpar(col = "darkgrey", fontsize = 10),
                                                x = unit(.2, "lines"), 
                                                y = unit(5, "lines"))),
                              layout_matrix = rbind(
                                c(1,2)),
                              widths = c(40,1))
    
    ggsave(file = paste0(png_chr_i_prefix, i, ".png"),
           plot = p_chr_i_arrow,
           height = 5, 
           width = 7.5, 
           units = "in",
           device = "png")
    
    if (i == 4 & zea == "maize") { #only chr 4 maize is a main figure, the others are supplemental
      ggsave(file = paste0(png_chr_i_prefix_main_lzw, i, ".tif"),
           plot = p_chr_i_arrow,
           height = 5, 
           width = 7.5, 
           units = "in",
           device = "tiff",
           dpi = 300, 
           compression = "lzw", type = "cairo")
    }else{
      ggsave(file = paste0(png_chr_i_prefix_supp_lzw, i, ".tif"),
             plot = p_chr_i_arrow,
             height = 5, 
             width = 7.5, 
             units = "in",
             device = "tiff",
             dpi = 300, 
             compression = "lzw", type = "cairo")
    }
  }
}



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


# ---------plot network of all observed data only -----------------
p_net_maize_data <- ggraph(peak_network[["maize"]], layout = "linear") +
  geom_edge_arc(aes(width = p_snps_shared_data*100,
                    alpha = p_snps_shared_data*100)) +
  geom_node_point(aes(color = ELEVATION), size = 3) +
  theme_graph(base_family = 'Helvetica') +
  geom_node_text(aes(label = name), 
                 angle = 90,
                 hjust = 0.5, #1,
                 y = -1.8, #-0.5,
                 repel = F) +
  scale_edge_alpha(range = c(0, 1), limits = c(0, 3.5)) +
  scale_edge_width(range = c(0, 3), limits = c(0, 3.5)) +
  scale_color_viridis(direction = -1) +
  coord_cartesian(clip = "off") +
  labs(edge_width = "% SNPs in\n shared peaks",
       edge_alpha = "% SNPs in\n shared peaks",
       color = "Elevation (m)") +
  ggtitle("Shared introgression peaks in maize")
# p_net_maize_data

p_net_mexicana_data <- ggraph(peak_network[["mexicana"]], layout = "linear") +
  geom_edge_arc(aes(width = p_snps_shared_data*100, 
                    alpha = p_snps_shared_data*100)) +
  geom_node_point(aes(color = ELEVATION), size = 3) + 
  theme_graph(base_family = 'Helvetica') +
  scale_edge_alpha(range = c(0, 1), limits = c(0, 3.5)) +
  scale_edge_width(range = c(0, 3), limits = c(0, 3.5)) +
  scale_color_viridis(direction = -1) +
  coord_cartesian(clip = "off") +
  labs(edge_width = "% SNPs in\n shared peaks",
       edge_alpha = "% SNPs in\n shared peaks",
       color = "Elevation (m)") +
  ggtitle("Shared introgression peaks in mexicana")
# p_net_mexicana_data

# combine maize and mexicana networks into 1 plot:
p_net_multi_data <- grid.arrange(grobs = list(ggplotGrob(p_net_maize_data +
                                                           theme(plot.margin = margin(c(t = 0, r = 5, b = 95, l = 2.5), unit = "pt"),
                                                                 plot.title = element_blank(),
                                                                 legend.position = "none"
                                                           )),
                                              ggplotGrob(p_net_mexicana_data +
                                                           theme(plot.margin = margin(c(t = 5, r = 5, b = 5, l = 2.5), unit = "pt"),
                                                                 plot.title = element_blank(),
                                                                 legend.position = "none"
                                                           )),
                                              cowplot::get_legend(p_net_maize_data),
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

# p_net_multi_data

ggsave(file = png_net_multi_data,
       plot = p_net_multi_data,
       height = 6.5, width = 6.2, 
       units = "in",
       device = "png") 
ggsave(file = png_net_multi_data_lzw,
       plot = p_net_multi_data,
       height = 6.5, width = 6.2, 
       dpi = 300, units = "in",
       device = "tiff",
       compression = "lzw", type = "cairo")


# ----------which pops commonly share peaks? plot ggupset or other combination matrix ---------
maize_pops_shared <- anc_outliers_list[["maize"]] %>% 
  arrange(ELEVATION) %>%
  mutate(LOCALITY_ELEVATION = paste0(LOCALITY,
                                     " ",
                                     ELEVATION,
                                     "m"),
         LOCALITY = reorder(LOCALITY, ELEVATION)) %>%
  filter(top_sd2) %>% 
  group_by(chr, pos) %>% 
  summarize(populations = paste0(sort(LOCALITY), collapse = "-"),
            n = length(unique(LOCALITY)))

p_combmatrix_maize <- maize_pops_shared %>%
  group_by(populations, n) %>%
  summarise(freq = n()/nrow(sites)*100) %>%
  ungroup() %>%
  mutate(populations = reorder(populations, -freq)) %>%
  arrange(populations) %>%
  head(n = 75) %>%
  ggplot(data = ., mapping = aes(x = populations, y = freq, fill = n)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  ggupset::axis_combmatrix(sep = "-", 
                  levels = meta_pops_list[["maize"]] %>%
                    arrange(desc(ELEVATION)) %>%
                    .$LOCALITY) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme_combmatrix(combmatrix.panel.point.color.fill = "#00BFC4",
                   combmatrix.panel.point.color.empty = "grey92",
                   combmatrix.panel.line.size = 0,
                   combmatrix.panel.point.size = 1.5,
                   combmatrix.panel.striped_background.color.one = "white",
                   combmatrix.panel.striped_background.color.two = "white",
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
ggsave(file = png_combmatrix_maize_lzw,
       plot = p_combmatrix_maize,
       height = 5, width = 7.5, 
       device = "tiff",
       dpi = 300, units = "in",
       compression = "lzw", type = "cairo")

# now mexicana
mexicana_pops_shared <- anc_outliers_list[["mexicana"]] %>% 
  arrange(ELEVATION) %>%
  mutate(LOCALITY_ELEVATION = paste0(LOCALITY,
                                     " ",
                                     ELEVATION,
                                     "m"),
         LOCALITY = reorder(LOCALITY, ELEVATION)) %>%
  filter(top_sd2) %>% 
  group_by(chr, pos) %>% 
  summarize(populations = paste0(sort(LOCALITY), collapse = "-"),
            n = length(unique(LOCALITY)))

p_combmatrix_mexicana <- mexicana_pops_shared %>%
  group_by(populations, n) %>%
  summarise(freq = n()/nrow(sites)*100) %>%
  ungroup() %>%
  mutate(populations = reorder(populations, -freq)) %>%
  arrange(populations) %>%
  head(n = 75) %>%
  ggplot(data = ., mapping = aes(x = populations, y = freq, fill = n)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  ggupset::axis_combmatrix(sep = "-", 
                  levels = meta_pops_list[["mexicana"]] %>%
                    arrange(desc(ELEVATION)) %>%
                    .$LOCALITY) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme_combmatrix(combmatrix.panel.point.color.fill = "#00BFC4",
                   combmatrix.panel.point.color.empty = "grey92",
                   combmatrix.panel.line.size = 0,
                   combmatrix.panel.point.size = 1.5,
                   combmatrix.panel.striped_background.color.one = "white",
                   combmatrix.panel.striped_background.color.two = "white",
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
ggsave(file = png_combmatrix_mexicana_lzw,
       plot = p_combmatrix_mexicana,
       height = 5, width = 7.5, 
       units = "in",
       device = "tiff",
       dpi = 300,
       compression = "lzw", type = "cairo") 


# ----------plot where shared peaks are across the genome-------
# introgression peaks defined by 2 sd > mean 
peaks_on_genome <- left_join(sites, 
                             bind_rows(outliers_by_pop_list[["maize"]],
                                       outliers_by_pop_list[["mexicana"]]),
          by = c("chr", "pos", "inv4m"))

top_shared_pops <- maize_pops_shared %>%
  group_by(populations, n) %>%
  summarise(freq = n()/nrow(sites)*100) %>%
  ungroup() %>%
  mutate(populations = reorder(populations, -freq)) %>%
  arrange(populations) %>%
  head(n = 100) %>%
  dplyr::select(populations, freq) 

p_peaks_on_genome <- top_shared_pops %>%
  left_join(., maize_pops_shared, by = "populations") %>%
  filter(n >= 5) %>%
  mutate(populations = reorder(populations, freq)) %>%
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

ggsave(file = png_peaks_on_genome,
       plot = p_peaks_on_genome,
       height = 5, width = 7.5, 
       units = "in",
       device = "png") 

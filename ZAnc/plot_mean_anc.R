#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)

# this script plots mean ancestry across the genome

# load variables from Snakefile
# prefix = "HILO_MAIZE55_PARV50"
# K = 3
zea = snakemake@params[["zea"]]
# zea = "maize"
ancestry = snakemake@params[["ancestry"]]
# ancestry = "mexicana"
colors_file = snakemake@input[["colors"]]
# colors_file = "colors.R"
anc_file = snakemake@input[["anc"]]
# anc_file = paste0("local_ancestry/results/ancestry_hmm/", prefix, "/K", K, "/Ne10000_yesBoot/anc/", zea, ".pops.anc.RData")
meta_file = snakemake@input[["meta_pop"]]
# meta_file = paste0("local_ancestry/results/ancestry_hmm/", prefix, "/K", K, "/Ne10000_yesBoot/anc/", zea, ".pop.meta.RData")
sites_file = snakemake@input[["sites"]]
# sites_file = paste0("local_ancestry/results/thinnedSNPs/", prefix, "/K", K, "/whole_genome.var.sites")
genome_file = snakemake@input[["genome"]]
# genome_file = "data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"
centromeres_file = snakemake@input[["centromeres"]]
# centromeres_file = "data/refMaize/centromere_positions_v4.txt"
fdr_in = snakemake@input[["fdr"]]
# fdr_in = paste0("ZAnc/results/", prefix, "/K", K, "/Ne10000_yesBoot/", zea, ".meanAnc.fdr.RData")
png_out = snakemake@output[["png"]]
# png_out = paste0("ZAnc/plots/", prefix, "_K", K, "_Ne10000_yesBoot_", zea, "_mean_", ancestry, "_anc.png")
rds = snakemake@output[["rds"]]
# rds = paste0("ZAnc/plots/", prefix, "_K", K, "_Ne10000_yesBoot_", zea, "_mean_", ancestry, "_anc.plot.rds")


# load data
source(colors_file)
load(anc_file)
load(meta_file)
load(fdr_in) # false discovery rates, FDRs

# load centromere positons
centromeres <- read.table(centromeres_file, header = T, stringsAsFactors = F,
                          sep = "\t")

# load chromosome lengths
genome <- read.table(genome_file, header = F, stringsAsFactors = F,
                     sep = "\t") %>%
  data.table::setnames(c("chr", "length")) %>%
  dplyr::mutate(chr_end = cumsum(length),
                chr_start = c(0, chr_end[1:(nrow(.)-1)])) %>%
  left_join(.,
            centromeres %>%
              dplyr::group_by(chr) %>% # chr9 has 2 segments that map to centromere region, so I get an approximate midpoint using both pieces
              summarise(cent_mid = 10^6 * (max(end) + min(start))/2), # convert from Mb to bp
            by = "chr") %>%
  dplyr::mutate(centromere = chr_start + cent_mid) # cumulative centromere positions

# and site/position information for SNPs with ancestry calls
sites <- read.table(sites_file, header = F, stringsAsFactors = F,
                    sep = "\t") %>%
  data.table::setnames(c("chr", "pos", "major", "minor")) %>%
  left_join(., genome, by = "chr") %>%
  dplyr::mutate(pos_cum = chr_start + pos) # get cumulative chromosomal position

# outlier plot whole genome
p_combined = bind_cols(sites, anc = anc_mean[[ancestry]]) %>%
  mutate(even_chr = ifelse(chr %% 2 == 0, "even", "odd"),
         zea = zea) %>%
  ggplot(.) +
  geom_hline(data = data.frame(intercept = filter(FDRs[[ancestry]], FDR == 0.05, threshold < Inf &
                                   threshold > -Inf)$threshold,
                                label = "fdr5",
                                stringsAsFactors = F),
              linetype = "solid",
              aes(yintercept = intercept,
                  color = label)) +
  geom_point(size = .1,
             aes(pos_cum, anc,
                 color = even_chr)) +
  geom_hline(data = data.frame(genomewide_mean = mean(anc_mean[[ancestry]]),
                                label = "genomewide_mean",
                                stringsAsFactors = F),
                aes(yintercept = genomewide_mean,
                    color = label),
              linetype = "dashed") +
  xlab("bp position on chromosomes (total length = 2.3Gb)") +
  ylab(paste("mean", ancestry, "ancestry")) +
  scale_colour_manual(values = c(odd = "darkgrey",
                                 even = unname(col_maize_mex_parv[zea]),
                                 genomewide_mean = "black",
                                 fdr5 = "#00BFC4"),
                      labels = c("odd chr", "even chr", "genomewide mean", "5% FDR"),
                      limits = c("odd", "even", "genomewide_mean", "fdr5")) +
  scale_x_continuous(label = genome$chr,
                     breaks = genome$centromere,
                     expand = expansion(mult = c(0, 0),
                                        add = c(0, 0))) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  ggtitle(paste("Sympatric", zea)) +
  ylim(0:1) +
  guides(color = guide_legend(override.aes = list(shape = NA,
                                                  linetype = c("dotted", "dotted", "dashed", "solid"),
                                                  width = 3)))
# p_combined

# save png
ggsave(plot = p_combined,
       file = png_out,
       height = 3, width = 12,
       units = "in", dpi = 300,
       device = "png")

# also save ggplot as R object
saveRDS(object = p_combined, file = rds)

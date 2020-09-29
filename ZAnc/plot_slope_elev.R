#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)

# this script plots slope from lm ancestry ~ elev across the genome,
# and calculates FDRs to show outlier loci

# load variables from Snakefile
zea = snakemake@params[["zea"]]
# zea = "maize"
colors_file = snakemake@input[["colors"]]
# colors_file = "colors.R"
fdr_file = snakemake@input[["fdr"]]
# fdr_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".lmElev.fdr.RData")
fit_file = snakemake@input[["fit"]]
# fit_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".lmElev.fit.RData")
meta_file = snakemake@input[["meta_pop"]]
# meta_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pop.meta.RData")
sites_file = snakemake@input[["sites"]]
# sites_file = "local_ancestry/results/thinnedSNPs/HILO_MAIZE55/whole_genome.var.sites"
genome_file = snakemake@input[["genome"]]
# genome_file = "data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"
png_out = snakemake@output[["png"]]
# png_out = paste0("ZAnc/plots/Ne10000_yesBoot/", zea, "_slope_elev.png")

# load data
source(colors_file)
load(fdr_file)
load(fit_file)
load(meta_file)

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
  left_join(., genome, by = "chr") %>%
  dplyr::mutate(pos_cum = chr_start + pos) # get cumulative chromosomal position

# mean of each chromosome (cumulative position)
axis_spacing = sites %>%
  group_by(chr) %>% 
  summarize(center=(max(pos_cum) + min(pos_cum)) / 2,
            start = min(pos_cum),
            end = max(pos_cum))

# outlier plot whole genome
p_elev = bind_cols(sites, fits) %>%
  mutate(even_chr = ifelse(chr %% 2 == 0, "even", "odd"),
         zea = zea) %>%
  ggplot(., aes(pos_cum, envWeights, 
                color = even_chr)) +
  geom_abline(slope = 0, intercept = filter(FDRs, FDR == 0.05)$thesholds, linetype = "solid", color = "blue") +
  geom_point(size = .1) +
  geom_abline(slope = 0, intercept = mean(fits$envWeights), color = "black", linetype = "dashed") +
  xlab("bp position on chromosomes (total length = 2.3Gb)") +
  ylab("slope ancestry ~ elev") +
  scale_colour_manual(values = c(odd = "darkgrey", even = unname(col_maize_mex_parv[zea]))) + 
  scale_x_continuous(label = axis_spacing$chr, breaks= axis_spacing$center) +
  theme(legend.position = "none") +
  theme_classic() +
  ggtitle(paste("Sympatric", zea, "- Change in mexicana ancestry over 1 km elevation gain")) +
  guides(color = F)
# p_elev

ggsave(plot = p_elev,
       file = png_out,
       height = 3, width = 12, 
       units = "in", dpi = 300,
       device = "png")

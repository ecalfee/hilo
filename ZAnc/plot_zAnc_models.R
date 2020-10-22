#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)

# this script plots results of all pop selection zAnc model
# and zElev model
# with FDRs to show outlier loci

# load variables from Snakefile
zea = snakemake@params[["zea"]]
# zea = "maize"
fdr_file = snakemake@input[["fdr"]]
# fdr_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".zAnc.fdr.RData")
colors_file = snakemake@input[["colors"]]
# colors_file = "colors.R"
fit_file = snakemake@input[["fit"]]
# fit_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".zAnc.fit.RData")
meta_file = snakemake@input[["meta_pop"]]
# meta_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pop.meta.RData")
sites_file = snakemake@input[["sites"]]
# sites_file = "local_ancestry/results/thinnedSNPs/HILO_MAIZE55/whole_genome.var.sites"
genome_file = snakemake@input[["genome"]]
# genome_file = "data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"

# png_out = snakemake@output[["png"]]
# png_out = paste0("ZAnc/plots/Ne10000_yesBoot/", zea, "_mean_anc.png")

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

# plot outliers for high mexicana ancestry
bind_cols(sites, fits) %>%
  filter(diff_all_ll2 > FDRs$thesholds[FDRs$tail == "high" & FDRs$FDR == 0.05 & FDRs$model == "all"]) %>%
  mutate(even_chr = ifelse(chr %% 2 == 0, "even", "odd"),
         zea = zea) %>%
  ggplot(., aes(pos_cum, diff_all_ll2, 
                color = even_chr)) +
  geom_point(size = .1) +
  xlim(c(0, max(sites$pos_cum)))
# and slope with elevation
bind_cols(sites, fits) %>%
  filter(diff_elev_ll2 > FDRs$thesholds[FDRs$tail == "high" & FDRs$FDR == 0.05 & FDRs$model == "elev"]) %>%
  mutate(even_chr = ifelse(chr %% 2 == 0, "even", "odd"),
         zea = zea) %>%
  ggplot(., aes(pos_cum, b1.elev, 
                color = even_chr)) +
  geom_point(size = .1) +
  xlim(c(0, max(sites$pos_cum)))


# outlier plot whole genome
p_combined = bind_cols(sites, anc = anc_mean) %>%
  mutate(even_chr = ifelse(chr %% 2 == 0, "even", "odd"),
         zea = zea) %>%
  ggplot(., aes(pos_cum, anc, 
                color = even_chr)) +
  geom_abline(slope = 0, intercept = filter(FDRs, FDR == 0.05, thesholds < Inf & thesholds > -Inf)$thesholds, linetype = "solid", color = "blue") +
  geom_point(size = .1) +
  geom_abline(slope = 0, intercept = mean(anc_mean), color = "black", linetype = "dashed") +
  xlab("bp position on chromosomes (total length = 2.3Gb)") +
  ylab("mean mexicana ancestry") +
  scale_colour_manual(values = c(odd = "darkgrey", even = unname(col_maize_mex_parv[zea]))) + 
  scale_x_continuous(label = axis_spacing$chr, breaks= axis_spacing$center) +
  theme(legend.position = "none") +
  theme_classic() +
  ggtitle(paste("Sympatric", zea)) +
  guides(color = F) +
  ylim(0:1)
# p_combined

ggsave(plot = p_combined,
       file = png_out,
       height = 3, width = 12, 
       units = "in", dpi = 300,
       device = "png")

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
meta_pops_file = snakemake@input[["meta_pop"]]
# meta_pops_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pop.meta.RData")
sites_file = snakemake@input[["sites"]]
# sites_file = "local_ancestry/results/thinnedSNPs/HILO_MAIZE55/whole_genome.var.sites"
png_out = snakemake@output[["png"]]
# png_out = paste0("mhl1_inv/plots/HILO_MAIZE55/Ne10000_yesBoot/mhl1_inv.png")
mhl1_bed = snakemake@input[["mhl1_bed"]]
# mhl1_bed = "data/known_QTL/chr9_bin4_mhl1_locus_v4.bed"


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

qtl9 <- read.table(mhl1_bed,
                   header = F, sep = "\t")
colnames(qtl9) <- c("chr", "start", "end")
bind_cols(sites, fits) %>%
  filter(chr == 9) %>%
  mutate(even_chr = ifelse(chr %% 2 == 0, "even", "odd"),
         zea = zea) %>%
  ggplot(., aes(pos/10^6, envWeights, 
                color = even_chr)) +
  geom_vline(xintercept = c(qtl9[["start"]], qtl9[["end"]])/10^6, linetype = "dashed", color = "orange") +
  geom_hline(yintercept = filter(FDRs, FDR == 0.05)$thesholds, linetype = "solid", color = "#00BFC4") +
  geom_point(size = .1) +
  geom_hline(yintercept = mean(fits$envWeights), color = "black", linetype = "dashed") +
  xlab("position on chr9 (Mbp)") +
  ylab("slope mexicana anc ~ elev") +
  ylim(c(-0.7, 1.2)) +
  scale_colour_manual(values = c(odd = "darkgrey", 
                                 even = unname(col_maize_mex_parv[zea]))) + 
  scale_x_continuous(expand = expansion(mult = c(0, 0),
                                        add = c(0, 0))) +
  theme(legend.position = "none") +
  theme_classic() +
  guides(color = F)


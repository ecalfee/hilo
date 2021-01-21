#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)

# this script plots mean ancestry across the genome,
# and calculates FDRs to show outlier loci

# load variables from Snakefile
zea = snakemake@params[["zea"]]
# zea = "maize"
fdr_functions = snakemake@input[["fdr_functions"]]
# fdr_functions = "ZAnc/FDR.R"
colors_file = snakemake@input[["colors"]]
# colors_file = "colors.R"
sim_file = snakemake@input[["sim"]]
# sim_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".MVN.RData")
anc_file = snakemake@input[["anc"]]
# anc_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pops.anc.RData")
meta_file = snakemake@input[["meta_pop"]]
# meta_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pop.meta.RData")
sites_file = snakemake@input[["sites"]]
# sites_file = "local_ancestry/results/thinnedSNPs/HILO_MAIZE55/whole_genome.var.sites"
genome_file = snakemake@input[["genome"]]
# genome_file = "data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"
png_out = snakemake@output[["png"]]
# png_out = paste0("ZAnc/plots/Ne10000_yesBoot/", zea, "_mean_anc.png")
fdr_out = snakemake@output[["fdr"]]
# fdr_out = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".meanAnc.fdr.RData")
rds_out = snakemake@output[["rds"]]
# rds = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".meanAnc.plot.rds")

# load data
source(fdr_functions)
source(colors_file)
load(anc_file)
load(sim_file)
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

# calculate false discovery rate thresholds
FDRs = calc_FDR(d = anc_mean, 
                s = MVN_mean, 
                test_values = seq(0, 1, by = .0001))

# what % of SNPs exceed these thresholds?
FDRs$n_SNPs = sapply(1:nrow(FDRs), function(i) 
  ifelse(FDRs$tail[i] == "high", 
         sum(anc_mean > FDRs$thesholds[i]),
         sum(anc_mean < FDRs$thesholds[i])))
FDRs$prop_SNPs = FDRs$n_SNPs/length(anc_mean)

save(FDRs, file = fdr_out)

# outlier plot whole genome
p_combined = bind_cols(sites, anc = anc_mean) %>%
  mutate(even_chr = ifelse(chr %% 2 == 0, "even", "odd"),
         zea = zea) %>%
  ggplot(.) +
  geom_hline(data = data.frame(intercept = filter(FDRs, FDR == 0.05, thesholds < Inf & 
                                   thesholds > -Inf)$thesholds, 
                                label = "fdr5",
                                stringsAsFactors = F),
              linetype = "solid", 
              aes(yintercept = intercept, 
                  color = label)) +
  geom_point(size = .1,
             aes(pos_cum, anc, 
                 color = even_chr)) +
  geom_hline(data = data.frame(genomewide_mean = mean(anc_mean),
                                label = "genomewide_mean",
                                stringsAsFactors = F),
                aes(yintercept = genomewide_mean, 
                    color = label), 
              linetype = "dashed") +
  xlab("bp position on chromosomes (total length = 2.3Gb)") +
  ylab("mean mexicana ancestry") +
  scale_colour_manual(values = c(odd = "darkgrey", 
                                 even = unname(col_maize_mex_parv[zea]),
                                 genomewide_mean = "black",
                                 fdr5 = "#00BFC4"),
                      labels = c("odd chr", "even chr", "genomewide mean", "5% FDR"),
                      limits = c("odd", "even", "genomewide_mean", "fdr5")) + 
  scale_x_continuous(label = axis_spacing$chr, 
                     breaks = axis_spacing$center,
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

ggsave(plot = p_combined,
       file = png_out,
       height = 3, width = 12, 
       units = "in", dpi = 300,
       device = "png")

# also save ggplot as R object
saveRDS(object = p_combined, file = rds)

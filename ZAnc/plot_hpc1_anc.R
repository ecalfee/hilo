#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(tidyr)

# this script plots slope from lm ancestry ~ elev
# and mexicana introgression into sympatric maize
# across the flowering time pathway gene HPC1

# load variables from Snakefile
zea = "maize"
colors_file = snakemake@input[["colors"]]
# colors_file = "colors.R"
fdr_lmElev_file = snakemake@input[["fdr_lmElev"]]
# fdr_lmElev_file = paste0("ZAnc/results/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/", zea, ".lmElev.fdr.RData")
fit_lmElev_file = snakemake@input[["fit_lmElev"]]
# fit_lmElev_file = paste0("ZAnc/results/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/", zea, ".lmElev.fit.RData")
fdr_maize_file = snakemake@input[["fdr_maize"]]
# fdr_maize_file = paste0("ZAnc/results/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/", zea, ".meanAnc.fdr.RData")
anc_maize_file = snakemake@input[["anc_maize"]]
# anc_maize_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/anc/", zea, ".pops.anc.RData")
meta_maize_file = snakemake@input[["meta_maize"]]
# meta_maize_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/anc/", zea, ".pop.meta.RData")
sites_file = snakemake@input[["sites"]]
# sites_file = "local_ancestry/results/thinnedSNPs/HILO_MAIZE55_PARV50/K3/whole_genome.var.sites"
png_out = snakemake@output[["png"]]
# png_out = "ZAnc/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_hpc1_lm_ancestry.png"
png_out_lzw = snakemake@output[["png_lzw"]]
# png_out_lzw = "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_hpc1_lm_ancestry.tif"
png_out_mex_anc = snakemake@output[["png_mex_anc"]]
# png_out_mex_anc = "ZAnc/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_hpc1_mean_ancestry.png"

# load data
source(colors_file)
load(fit_lmElev_file)
load(fdr_lmElev_file)
assign("FDRs_lmElev", FDRs)
load(meta_maize_file)
load(anc_maize_file)
load(fdr_maize_file)
assign("FDRs_meanAnc", FDRs)
FDRs <- bind_rows(mutate(FDRs_lmElev, model = "lmElev"),
                  mutate(FDRs_meanAnc[["mexicana"]], model = "meanAnc")) %>%
  mutate(ancestry = "mexicana")


# load site/position information for SNPs with ancestry calls
sites <- read.table(sites_file, header = F, stringsAsFactors = F,
                    sep = "\t") %>%
  data.table::setnames(c("chr", "pos", "major", "minor"))

hpc1 <- list(gene = "hpc1", chr = 3,
                   start = 7735746,
                   end = 7737263)

# plot mean mexicana introgression around hpc1
p_hpc1_anc <- mutate(sites, mexicana_anc = anc_mean[["mexicana"]]) %>%
  filter(chr == 3, pos/10^6 > 5, pos/10^6 < 11) %>%
  ggplot(., aes(pos/10^6, mexicana_anc)) +
  geom_hline(yintercept = filter(FDRs, model == "meanAnc", ancestry == "mexicana")$threshold, linetype = "solid", color = "#00BFC4") +
  geom_point(size = .1, color = "darkgrey") +
  geom_hline(yintercept = mean(fits$envWeights), color = "black", linetype = "dashed") +
  xlab("position on chr3 (Mbp)") +
  ylab("slope mexicana ancestry ~ elevation") +
  scale_y_continuous(breaks = c(-.5, 0, .5, 1),
                     limits = c(-0.7, 1.45),
                     expand = expansion(mult = c(0, 0),
                                        add = c(0, 0))) +
  theme(legend.position = "none") +
  theme_classic() +
  geom_text(label = "HPC1", x = sum(c(hpc1[["start"]], hpc1[["end"]]))/(2*10^6), size = 4,
            color = "black", y = 1.37) +
  geom_segment(x = hpc1[["start"]]/10^6, xend = hpc1[["end"]]/10^6, y = 1.21, yend = 1.21, color = "black") +
  geom_segment(x = hpc1[["start"]]/10^6, xend = hpc1[["start"]]/10^6, y = 1.17, yend = 1.25, color = "black") +
  geom_segment(x = hpc1[["end"]]/10^6, xend = hpc1[["end"]]/10^6, y = 1.17, yend = 1.25, color = "black")
#p_hpc1_anc

# make simple plot just to extract legend only for genomewide mean and 5% FDR lines
plot_for_legend_only <- ggplot(data = data.frame(lines = c("genomewide_mean", "fdr5"), 
                                                 value = 1)) +
  geom_hline(aes(yintercept = value, color = lines, linetype = lines)) +
  theme_classic() +
  labs(color = "a", linetype = "a") +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  scale_colour_manual(values = c(genomewide_mean = "black",
                                 fdr5 = "#00BFC4"),
                      labels = c("genomewide mean", "10%, 5%, 1% FDR"),
                      limits = c("genomewide_mean", "fdr5")) +
  scale_linetype_manual(values = c(genomewide_mean = "dashed",
                                   fdr5 = "solid"),
                        labels = c("genomewide mean", "10%, 5%, 1% FDR"),
                        limits = c("genomewide_mean", "fdr5"))

p_hpc1_anc_legend <- grid.arrange(grobs = list(ggplotGrob(p_hpc1_anc),
                               cowplot::get_legend(plot_for_legend_only)),
                  nrow = 2,
                  heights = c(1,.1))
ggsave(png_out_mex_anc, 
       plot = p_hpc1_anc_legend, 
       device = "png", 
       width = 7.5, 
       height = 4, 
       units = "in",
       dpi = 300)

# plot slope across elevation around hpc1
p_hpc1_lm <- bind_cols(sites, fits) %>%
  filter(chr == 3, pos/10^6 > 5, pos/10^6 < 11) %>%
  ggplot(., aes(pos/10^6, envWeights)) +
  geom_hline(yintercept = filter(FDRs, model == "lmElev", tail == "high")$threshold, linetype = "solid", color = "#00BFC4") +
  geom_point(size = .1, color = "darkgrey") +
  geom_hline(yintercept = mean(fits$envWeights), color = "black", linetype = "dashed") +
  xlab("position on chr3 (Mbp)") +
  ylab("slope mexicana ancestry ~ elevation") +
  scale_y_continuous(breaks = c(-.5, 0, .5, 1),
                     limits = c(-0.7, 1.45),
                     expand = expansion(mult = c(0, 0),
                                        add = c(0, 0))) +
  theme(legend.position = "none") +
  theme_classic() +
  geom_text(label = "HPC1", x = sum(c(hpc1[["start"]], hpc1[["end"]]))/(2*10^6), size = 4,
            color = "black", y = 1.37) +
  geom_segment(x = hpc1[["start"]]/10^6, xend = hpc1[["end"]]/10^6, y = 1.21, yend = 1.21, color = "black") +
  geom_segment(x = hpc1[["start"]]/10^6, xend = hpc1[["start"]]/10^6, y = 1.17, yend = 1.25, color = "black") +
  geom_segment(x = hpc1[["end"]]/10^6, xend = hpc1[["end"]]/10^6, y = 1.17, yend = 1.25, color = "black")
#p_hpc1_lm
# at a few snps with the highest slopes mexicana anc ~ elevation, plot ancestry for each population across elevation
top_snps <- bind_cols(sites, fits) %>%
  filter(chr == 3, pos/10^6 > 5, pos/10^6 < 12) %>%
  mutate(hpc1_20kb = (pos + 20000 > hpc1[["start"]] & pos - 20000 < hpc1[["end"]]),
         in_hpc1 = (pos > hpc1[["start"]] & pos < hpc1[["end"]])) %>%
  filter(hpc1_20kb) %>%
  dplyr::arrange(., desc(envWeights)) %>%
  mutate(snp = 1:nrow(.)) %>%
  filter(snp <= 3)
# snp1 within 20kb of upstream of gene: 3 7733803     G     C   -1.204678  0.7961550
# snp2 within gene: 3 7737055     C     T   -1.183194  0.7871671
# snp3 within 20kb of downstream of gene: 3 7739077     A     G   -1.174751  0.7846878
a <- bind_cols(sites, data.frame(anc[["mexicana"]])) %>%
  left_join(top_snps, ., by = c("chr", "pos", "major", "minor"))
pop_means = apply(anc[["mexicana"]], 2, mean) %>%
  t(.) %>%
  data.frame(.) %>%
  pivot_longer(cols = starts_with("pop"), names_to = "pop", values_to = "mean_mexicana_anc") %>%
  dplyr::mutate(SNP = "genomewide mean") %>%
  left_join(., meta_pops, by = "pop")
p_snps = a %>%
  filter(in_hpc1) %>% # only plot top snp within hpc1
  pivot_longer(cols = starts_with("pop"), names_to = "pop", values_to = "mexicana_anc") %>%
  left_join(., meta_pops, by = "pop") %>%
  dplyr::mutate(SNP = paste(chr, pos, sep = ":")) %>%
  ggplot(., aes(x = ELEVATION, y = mexicana_anc, color = SNP)) +
  geom_point() +
  geom_abline(aes(intercept = `(Intercept)`, slope = envWeights/1000, color = paste(chr, pos, sep = ":"))) +
  geom_point(data = pop_means, aes(y = mean_mexicana_anc)) +
  theme_light() +
  ylab("Introgressed mexicana ancestry") +
  xlab("Elevation (km)") +
  labs(color = "SNP") +
  scale_color_manual(values = c(#"#D55E00", 
                                "#009E73", 
                                #"#F0E442", 
                                "Black"))


p_lm <- grid.arrange(grobs = list(textGrob(label = "A",
                                           just = "top",
                                           x = unit(0.5, "lines")),
                                  ggplotGrob(p_hpc1_lm),
                                  cowplot::get_legend(plot_for_legend_only),
                                  textGrob(label = "B",
                                           just = "top",
                                           x = unit(0.5, "lines")),
                                  ggplotGrob(p_snps)),
                  nrow = 5,
                  heights = c(.1,1,.1,.1,1))

# p_lm
ggsave(png_out, 
       plot = p_lm, 
       device = "png", 
       width = 7.5, 
       height = 7.5, 
       units = "in",
       dpi = 300)

ggsave(png_out_lzw, 
       plot = p_lm, 
       device = "tiff", 
       width = 7.5, 
       height = 7.5, 
       units = "in",
       dpi = 300,
       compression = "lzw", type = "cairo")
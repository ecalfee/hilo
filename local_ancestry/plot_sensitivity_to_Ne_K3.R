#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(GGally)
# this script plots timing of admixture estimated from ancestry_hmm for different Nes
# and correlations between local ancestry estimates for different Nes

options(scipen = 999) # don't convert to scientific notation

# load variables from Snakefile
prefix = snakemake@params[["prefix"]]
# prefix = "HILO_MAIZE55_PARV50"

png_times = snakemake@output[["png_times"]]
# png_times = paste0("local_ancestry/plots/", prefix, "_K3_sensitivity_to_Ne_admix_times.png")
png_times_lzw = snakemake@output[["png_times_lzw"]]
# png_times_lzw = paste0("../hilo_manuscript/figures_supp/", prefix, "_K3_sensitivity_to_Ne_admix_times.tif")
png_anc = snakemake@output[["png_anc"]]
# png_anc = paste0("local_ancestry/plots/", prefix, "_K3_sensitivity_to_Ne_local_ancestry.png")
png_anc_lzw = snakemake@output[["png_anc_lzw"]]
# png_anc_lzw = paste0("../hilo_manuscript/figures_supp/", prefix, "_K3_sensitivity_to_Ne_local_ancestry.tif")

# variables
Nes = c(1000, 10000, 100000)
zea = c("maize", "mexicana")
ancestries = c("maize", "mexicana", "parv")
YESNO = "no"

# load data
source(snakemake@params[["colors"]]) # plotting colors
# source("colors.R")
load(paste0("samples/", prefix, "_meta.RData"))

# plot timing estimates
times = list()

for (z in zea){
  # load meta data
  load(paste0("local_ancestry/results/ancestry_hmm/",
              prefix, "/K3/Ne", Ne, "_", YESNO,
              "Boot/anc/", z, ".pop.meta.RData"))
  times[[z]] = do.call(bind_rows, lapply(meta_pops$pop, function(p)
    do.call(bind_rows, lapply(Nes, function(N)
    read.table(paste0("local_ancestry/results/ancestry_hmm/",
                      prefix, "/K3/Ne", N, "_", YESNO,
                      "Boot/", p, ".times"),
               header = T, stringsAsFactors = F, sep = " ") %>%
      dplyr::mutate(pop = p,
                    Ne = N))))) %>%
    left_join(., dplyr::select(meta_pops, pop, zea, LOCALITY, ELEVATION), by = "pop")
}

p_times <- do.call(bind_rows, times) %>%
  rename(time = t) %>%
  mutate(Ne = factor(Ne)) %>%
  mutate(zea = paste("sympatric", zea)) %>%
  mutate(admixture_pulse = ifelse(ancestry == "parv", "parviglumis", ancestry)) %>%
  ggplot(., aes(x = reorder(LOCALITY, ELEVATION),
                y = time,
                shape = Ne,
                col = admixture_pulse)) +
  geom_point() + # plot time estimate
  xlab("Sympatric population") +
  ylab("Admixture time (generations)") +
  theme_classic() +
  ylim(-100, 1650) +
  guides(color = guide_legend("Admixture pulse"),
         alpha = guide_legend("Pulse size")) +
  scale_color_manual(values = col_maize_mex_parv) +
  scale_shape_manual(values = c(19, 5, 3)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~zea)

#p_times
ggsave(file = png_times,
       plot = p_times,
       device = "png",
       width = 5, height = 4,
       units = "in", dpi = 300)

ggsave(file = png_times_lzw,
       plot = p_times,
       device = "tiff",
       width = 5, height = 4,
       compression = "lzw", type = "cairo",
       units = "in", dpi = 300)

# plot correlations for mean local ancestry across the genome
# pairs plot
# load ancestry data
a = "parv"
z = "maize"
d <- do.call(bind_rows, lapply(Nes, function(N)
  read.table(paste0("local_ancestry/results/ancestry_hmm/", PREFIX, "/K3/Ne", N, "_noBoot/anc/", z, ".", a, "_anc.bed"),
                header = T) %>%
  mutate(Ne = paste0("Ne=", N), 
         ancestry = a, 
         zea = paste("sympatric", z)))) %>%
  dplyr::select(zea, ancestry, Ne, anc_freq) %>%
  group_by(Ne) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(data = ., names_from = Ne, values_from = anc_freq) %>%
  dplyr::select(-row)
ggpairs(d, 
        aes(color = ancestry, shape = zea))                
ggpairs(filter(d, zea == "sympatric maize"), 
        aes(color = ancestry)) +
  scale_color_manual(values = col_maize_mex_parv)
# what? this shouldn't be a correlation of 1 ...
# but I could just plot/list the correlations instead of a complicated pairs plot
# and maybe plot just mexicana ancestry in maize?
ggplot(d, aes(color = ancestry, shape = zea, x = `Ne=1000`, y = `Ne=100000`)) +
  geom_point() +
  geom_smooth(method = "lm")
#+
 # scale_shape_manual(values = shape_maize_mex_parv) +
                
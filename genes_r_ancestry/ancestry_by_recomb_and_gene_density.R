# in this script I use recombination rate
# and gene density to predict ancestry patterns across the genome
# testing the hypothesis that minor ancestry is lowest in
# genomic regions with high gene density and low recombination
library(quantreg) # for regression on quantiles (percent ancestry)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)

# first, load data:
#source("ZAnc_statistic.R") 
#d = data.frame(pos, ZAnc_maize = maize[["ZAnc"]], ZAnc_mex = mex[["ZAnc"]], ZAnc_all = all[["ZAnc"]])
PREFIX="pass2_alloMAIZE"
LOCAL_ANC_SUBDIR="output_noBoot"
dir_in = paste0("../local_ancestry/results/ancestry_hmm/", PREFIX)
dir_anc = file.path(dir_in, LOCAL_ANC_SUBDIR, "anc")
dir_sites = paste0("../local_ancestry/results/thinnedSNPs/", PREFIX)

# metadata for each individual, including pop association
pop_elev <- read.table("../data/riplasm/gps_and_elevation_for_sample_sites.txt",
                       stringsAsFactors = F, header = T, sep = "\t") %>%
  dplyr::select(., popN, ELEVATION)
hilo <- read.table("../samples/hilo_meta.txt", stringsAsFactors = F, header = T, sep = "\t")
# admixture proportions per pop
meta.pops <- read.table(paste0("../global_ancestry/results/NGSAdmix/", PREFIX, "/globalAdmixtureIncludedByPopN.txt"),
                        stringsAsFactors = F, header = T) %>%
  left_join(., pop_elev, by = c("popN")) %>%
  left_join(., unique(select(hilo, c("popN", "zea", "symp_allo", "RI_ACCESSION", "GEOCTY", "LOCALITY"))), by = c("popN")) %>%
  mutate(pop = paste0("pop", popN))

meta.ind <- read.table(paste0("../global_ancestry/results/NGSAdmix/", PREFIX, "/globalAdmixtureByIncludedIndividual.txt"), 
                   stringsAsFactors = F, 
                   header = T, sep = "\t") %>%
  left_join(., hilo, by = c("ID", "popN")) %>%
  left_join(., pop_elev, by = c("popN"))%>%
  mutate(pop = paste0("pop", popN))

# does mean coverage correlated with estimated genomewide indidivual ancestry?
# we might expect individuals with more mexicana to map less well and have on avg. slightly lower coverage
ggplot(meta.ind, aes(x = est_coverage, y = alpha_mex, color = pop)) +
  geom_point(alpha = .8)

# get inversions:
inv = read.table("../data/refMaize/inversions/knownInv_v4_coord.txt",
                 stringsAsFactors = F, header = T)
colnames(inv) = c("ID", "chr", "start", "end", "length")
excl_inv_buffer = 1000000 # exclude 1 megabase around inversion
inv$excl_start = inv$start - excl_inv_buffer
inv$excl_end = inv$end + excl_inv_buffer
# which inversions are segregating and therefore should be excluded from some analyses?
# this list of inversions is checked visually for long ancestry blocks at potential inversion locations
inv$present = (inv$ID %in% c("inv4m"))

# color palette for locations (n=13 locations)
location_colors = c("gray10", "deepskyblue1",
                    brewer.pal(11, "Spectral"))
# color palette for mex and maize
blues = brewer.pal(n = 9, name = "Blues")[c(6,8)]
yellows = brewer.pal(n = 9, name = "YlOrBr")[c(4,6)]
colors_maize2mex = c(yellows, blues)
labels_maize2mex = c("allopatric_maize", "sympatric_maize", "sympatric_mexicana", "allopatric_mexicana")

# read in population ancestry input files from calc_genomewide_pop_anc_freq.R
pop_anc_list = lapply(meta.pops$popN, function(pop) read.table(paste0(dir_anc, "/pop", pop, ".anc.freq"), 
                                                          stringsAsFactors = F))

# combine population ancestry frequencies into a matrix
# where rows are populations and columns are snps
all_anc = do.call(cbind, pop_anc_list)
colnames(all_anc) <- meta.pops$pop

# I'll get position info from var.sites used for ancestry_hmm
# and then gene density from the 0.1cM windows around each SNP
sites <- do.call(rbind, 
                 lapply(1:10, function(i)
                   read.table(paste0(dir_sites, "/chr", i, ".var.sites"),
                              header = F, stringsAsFactors = F)))
colnames(sites) <- c("chr", "pos", "major", "minor")
wind <- read.table(paste0(dir_sites, "/windows0.1cM/gene_overlap.bed"),
                   header = F, stringsAsFactors = F)
colnames(wind) <- c("chr_wind", "start", "end", 
                    "width_cM", "coding_bp",
                    "n_CDS", "width_bp", # note number CDS is number contigous coding sequences (may be more than one gene (?))
                    "perc_coding")
pos <- cbind(sites, wind)
# calculate local recombination rate
# in cM/Mb
pos$rate <- pos$width_cM/(pos$width_bp/10^6)
hist(pos$rate)
summary(pos$rate)
# quantiles for defining 'low' and 'high' recomb.
low_quant = .2
high_quant = .8
quantile(pos$rate, c(low_quant, high_quant))

# group rates by quantiles for barplots:
pos$bin_rate = cut(pos$rate, # include full range of depth in cut
                 breaks = quantile(pos$rate, p = seq(0, 1, by = .1)),
                 right = T,
                 include.lowest = T)
pos$coding_density = pos$coding_bp/pos$width_cM
pos$bin_coding_density = cut(pos$coding_density,
                           breaks = quantile(pos$coding_density, p = seq(0, 1, by = .2)),
                           right = T,
                           include.lowest = T)
# population mean alpha is the mean of each population's mean (so each pop is weighted equally indep. of its sample size)
# whereas mean alpha is the mean mexicana ancestry at a locus across all individuals sampled
# populations are roughly sampled equally, so they are highly correlated, but not the same
# make a summary directory d
d <- pos # temporary until I get 0.1cM windows analysis done
maize_pops <- unique(meta.pops$pop[meta.pops$zea == "maize"])
mexicana_pops <- unique(meta.pops$pop[meta.pops$zea == "mexicana"])
maize_anc <- all_anc[ , maize_pops]
mexicana_anc <- all_anc[ , mexicana_pops]
d$pop_meanAlpha_maize = rowMeans(all_anc[ , maize_pops])
d$pop_meanAlpha_mex = rowMeans(all_anc[ , mexicana_pops])
N_per_pop_mex <- sapply(mexicana_pops, function(pop) unique(meta.pops$n[meta.pops$pop == pop]))
d$meanAlpha_mex = t((N_per_pop_mex %*% t(all_anc[ , mexicana_pops]))/sum(N_per_pop_mex))[ , 1]
N_per_pop_maize <- sapply(mexicana_pops, function(pop) unique(meta.pops$n[meta.pops$pop == pop]))
d$meanAlpha_maize = t((N_per_pop_maize %*% t(all_anc[ , maize_pops]))/sum(N_per_pop_maize))[ , 1]

# within an inversion?
d <- mutate(d, inv4m = (chr == inv$chr[inv$ID=="inv4m"] & 
                          pos >= inv$excl_start[inv$ID=="inv4m"] & 
                          pos <= inv$excl_end[inv$ID == "inv4m"])) 
d <- mutate(d, inv9d_or_e = chr == inv$chr[inv$ID=="inv9d_or_e"] &
              pos >= inv$excl_start[inv$ID == "inv9d_or_e"] &
              pos <= inv$excl_end[inv$ID == "inv9d_or_e"])
d <- mutate(d, inv9e_or_d = chr == inv$chr[inv$ID == "inv9e_or_d"] &
              pos >= inv$excl_start[inv$ID == "inv9e_or_d"] &
              pos <= inv$excl_end[inv$ID == "inv9e_or_d"])

d2 <- d %>%
  gather(., "group_alpha", "alpha", 
         c("meanAlpha_mex", "meanAlpha_maize"))

# exclude inversion area
d_split <- d %>%
  filter(., !inv_any) %>%
  gather(., "group_ZAnc", "ZAnc",
         c("ZAnc_mex", "ZAnc_maize", "ZAnc_all"))

# make bar plot
box_rate_ZAnc <- 
  ggplot(d_split, aes(y = ZAnc, x = bin_rate,
                      fill = group_ZAnc)) + 
  scale_fill_manual(values = c("grey", colors_maize2mex[c(1,4)])) +
  geom_boxplot() +
  facet_wrap(~group_ZAnc) +
  labs(x = "recombination rate (0.1cM winds)", 
       y = "ZAnc statistic by group",
       main = "Effect of local recomb. on ZAnc statistic")
plot(box_rate_ZAnc)
ggsave("../plots/box_rate_ZAnc_0.1cM_windows.png", 
       height = 15, width = 20, 
       units = "in", device = "png")

# make bar plot for mean ancestry mexicna
box_rate_alpha <- d2 %>%
  filter(., !inv4m) %>%
  ggplot(., aes(y = alpha, x = bin_rate,
                      fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_boxplot() +
  facet_wrap(~group_alpha) +
  labs(x = "recombination rate (0.1cM winds)", 
       y = "Mean mex. ancestry by group",
       main = "Effect of local recomb. on mean mex ancestry (alpha)")
plot(box_rate_alpha)
ggsave("plots/box_rate_alpha_0.1cM_windows_10rbins.png", 
       height = 15, width = 20, 
       units = "in", device = "png")
# same but as a violin plot to show distribution:
violin_rate_alpha <- 
  d2 %>%
  filter(., !inv4m) %>%
  ggplot(., aes(y = alpha, x = bin_rate,
                      fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_violin() +
  facet_wrap(~group_alpha) +
  labs(x = "recombination rate (0.1cM winds)", 
       y = "Mean mex. ancestry by group",
       main = "Effect of local recomb. on mean mex ancestry (alpha)")
plot(violin_rate_alpha)
ggsave("plots/violin_rate_alpha_0.1cM_windows.png",
       height = 15, width = 20,
       units = "in", device = "png")
# rather than binning, I can also view this information in a scatter plot and regression framework
# WARNING: Loess takes a really long time (!) because there are so many data points
# I'll start with maize
png(paste("../plots/smooth_scatter_alpha_by_recomb_rate_maize.png"), # saves plot as pin in ../plots/
    height = 4, width = 7, units = "in", res = 150)
with(d_split[d_split$group_alpha=="meanAlpha_maize",], smoothScatter(x = rate, y = alpha, main = "mean mexicana ancestry in maize predicted by recomb. rate"))
#with(d_split[d_split$group_alpha=="meanAlpha_maize",], lines(loess(alpha ~ rate), col = "red")) # too slow, plot a subset:
j <- sample(1:nrow(d_split[d_split$group_alpha=="meanAlpha_maize",]), 100000, replace = F)
with(d_split[d_split$group_alpha=="meanAlpha_maize",], lines(lowess(alpha[j] ~ rate[j]), col = "red"))
dev.off()

# it's helpful to take the log of recombination rate because r is so variable and we might expect diminishing returns
png(paste("../plots/smooth_scatter_alpha_by_log_recomb_rate_maize.png"), # saves plot as pin in ../plots/
    height = 4, width = 7, units = "in", res = 150)
with(d_split[d_split$group_alpha=="meanAlpha_maize",], 
     smoothScatter(x = log(rate), y = alpha, main = "mean mexicana ancestry in maize predicted by log recomb. rate"))
#with(d_split[d_split$group_alpha=="meanAlpha_maize",], lines(loess(alpha ~ log(rate)), col = "red")) # too slow, plot a subset:
with(d_split[d_split$group_alpha=="meanAlpha_maize",], lines(lowess(alpha[j] ~ log(rate[j])), col = "red"))
dev.off()

# now mexicana
png(paste("../plots/smooth_scatter_alpha_by_recomb_rate_mex.png"), # saves plot as pin in ../plots/
    height = 4, width = 7, units = "in", res = 150)
with(d_split[d_split$group_alpha=="meanAlpha_mex",], 
     smoothScatter(x = rate, y = alpha, main = "mean mexicana ancestry in mexicana predicted by recomb. rate"))
jmex <- sample(1:nrow(d_split[d_split$group_alpha=="meanAlpha_mex",]), 1000000, replace = F)
with(d_split[d_split$group_alpha=="meanAlpha_mex",], lines(lowess(alpha[jmex] ~ rate[jmex]), col = "red"))
dev.off()
# log transformed:
png(paste("../plots/smooth_scatter_alpha_by_recomb_log_rate_mex.png"), # saves plot as pin in ../plots/
    height = 4, width = 7, units = "in", res = 150)
with(d_split[d_split$group_alpha=="meanAlpha_mex",], 
     smoothScatter(x = log(rate), y = alpha, main = "mean mexicana ancestry in mexicana predicted by log recomb. rate"))
with(d_split[d_split$group_alpha=="meanAlpha_mex",], lines(lowess(alpha[jmex] ~ log(rate)[jmex]), col = "red"))
dev.off()

# using ggplot2 (subset because o/w takes a long time):
ggplot(d_split[sample(1:nrow(d_split), 100000, replace = F),], aes(y = alpha, x = rate,
                                                                   color = group_alpha)) + 
  scale_color_manual(values = colors_maize2mex[c(1,4)]) +
  geom_point(alpha = .5) +
  facet_wrap(~group_alpha) +
  geom_smooth(color = "black") +
  labs(x = "recombination rate (0.1cM winds)", 
       y = "Mean mex. ancestry (by group)all in",
       main = "Effect of local recomb. on mean mex ancestry (alpha)")
ggsave("../plots/scatter_rate_alpha_0.1cM_windows.png",
       height = 15, width = 20,
       units = "in", device = "png")
ggplot(d_split[sample(1:nrow(d_split), 100000, replace = F),], aes(y = alpha, x = log(rate),
                                                                   color = group_alpha)) + 
  scale_color_manual(values = colors_maize2mex[c(1,4)]) +
  geom_point(alpha = .5) +
  facet_wrap(~group_alpha) +
  geom_smooth(color = "black") +
  labs(x = "log recombination rate (0.1cM winds)", 
       y = "Mean mex. ancestry (by group)all in",
       main = "Effect of local recomb. on mean mex ancestry (alpha)")
ggsave("../plots/scatter_log_rate_alpha_0.1cM_windows.png",
       height = 15, width = 20,
       units = "in", device = "png")


# now for gene density: make bar plot for mean ancestry mexicna
box_coding_density_alpha <- 
  d2 %>%
  filter(., !inv4m) %>%
  ggplot(., aes(y = alpha, x = bin_coding_density,
                      fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_boxplot() +
  facet_wrap(~group_alpha) +
  labs(x = "coding density (bp/cM)", 
       y = "Mean mex. ancestry by group",
       main = "Effect of local coding density on mean mex ancestry (alpha)")
plot(box_coding_density_alpha)
ggsave("plots/box_coding_density_alpha_0.1cM_windows.png", 
       height = 15, width = 20, 
       units = "in", device = "png")

# now for gene density: make bar plot for mean ancestry mexicna
violin_coding_density_alpha <- 
  d2 %>%
  filter(., !inv4m) %>%
  ggplot(., aes(y = alpha, x = bin_coding_density,
                      fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_violin() +
  facet_wrap(~group_alpha) +
  labs(x = "coding density (bp/cM)", 
       y = "Mean mex. ancestry by group",
       main = "Effect of local coding density on mean mex ancestry (alpha)")
plot(violin_coding_density_alpha)
ggsave("plots/violin_coding_density_alpha_0.1cM_windows.png", 
       height = 15, width = 20, 
       units = "in", device = "png")
# I can also look at alpha ~ coding density using a smoothed scatter plot
# in maize:
png(paste("../plots/smooth_scatter_alpha_by_coding_density_maize.png"), # saves plot as pin in ../plots/
    height = 4, width = 7, units = "in", res = 150)
with(d_split[d_split$group_alpha=="meanAlpha_maize",], 
     smoothScatter(x = coding_density, y = alpha, main = "mean mexicana ancestry in maize predicted by coding density"))
with(d_split[d_split$group_alpha=="meanAlpha_maize",], 
     lines(lowess(alpha[j] ~ coding_density[j]), col = "red"))
dev.off()

# and in mexicana:
png(paste("../plots/smooth_scatter_alpha_by_coding_density_mex.png"), # saves plot as pin in ../plots/
    height = 4, width = 7, units = "in", res = 150)
with(d_split[d_split$group_alpha=="meanAlpha_mex",], smoothScatter(x = coding_density, y = alpha, main = "mean mexicana ancestry in mexicana predicted by coding density"))
with(d_split[d_split$group_alpha=="meanAlpha_mex",], lines(lowess(alpha[jmex] ~ coding_density[jmex]), col = "red"))
dev.off()

# plot ZAnc against gene density in 0.1cM windows around SNPs
d %>%
  #filter(., chr == "4") %>%
  ggplot(., aes(x = perc_coding, y = ZAnc_maize)) + 
  geom_point(cex= .1, alpha = .5, col = "black") +
  geom_smooth(col = "blue")


# plot recomb. rate against gene density in 0.1cM windows around SNPs
d %>%
  #filter(., chr == "4") %>%
  ggplot(., aes(x = perc_coding, y = rate)) + 
  geom_point(cex= .1, alpha = .5, color = "black") +
  geom_smooth(col = "blue") + 
  ggtitle("higher recombination in gene-dense .1cM windows")
ggsave("../plots/scatter_recomb_rate_by_coding_density_0.1cM_windows.png", 
       height = 15, width = 20, 
       units = "in", device = "png")
d %>%
  #filter(., chr == "4") %>%
  ggplot(., aes(y = perc_coding, x = log(rate))) + 
  geom_point(cex= .1, alpha = .5, col = "black") +
  geom_smooth(col = "blue")

# I can describe in a linear model the effects of total coding bp and non-coding bp in a .1cM window:
lm_maize <- d %>%
  mutate(noncoding_bp = width_bp - coding_bp) %>%
  lm(meanAlpha_maize ~ coding_bp + noncoding_bp, data = .)
lm_mex <- d %>%
  mutate(noncoding_bp = width_bp - coding_bp) %>%
  lm(meanAlpha_mex ~ coding_bp + noncoding_bp, data = .)
summary(lm_maize)
summary(lm_mex)
# linear regression isn't right because alpha is bounded by 0-1
# so I'm doing quantile regression below using the quantreg package:
# I don't quite understand how tau is used
# pfn method is recommended for computational efficiency
rq_maize <- d %>%
  mutate(coding_kb = coding_bp/1000) %>%
  mutate(noncoding_kb = (width_bp - coding_bp)/1000) %>%
  quantreg::rq(meanAlpha_maize ~ coding_kb + noncoding_kb, tau = .5, data = ., method = "pfn")
summary(rq_maize)
rq_mex <- d %>%
  mutate(coding_kb = coding_bp/1000) %>%
  mutate(noncoding_kb = (width_bp - coding_bp)/1000) %>%
  quantreg::rq(meanAlpha_mex ~ coding_kb + noncoding_kb, tau = .5, data = ., method = "pfn")
summary(rq_mex)
# I can't explain why coding_bp would increase mexicana in mexicana BUT non-coding bp would decrease mexicana in mexicana
# unless there are perhaps key inversions with high recomb. rates favoring mexicana ancestry which I have overlooked;
# or maybe less likely low recombination regions with universally favored maize ancestry.
# If I split populations out by elevation, I could also test whether populations at lower elevations
# select more strongly against mexicana with a coding_bp*elevation term

# these lm's suggest a better predictor should be the number of coding bp within .1cM
# plot this (takes a while to plot and to save):
d_split %>%
#d_split[sample(1:nrow(d_split), 100000, replace = F),] %>%
  ggplot(., 
       aes(y = alpha, x = coding_bp, 
           color = group_alpha)) + 
  scale_color_manual(values = colors_maize2mex[c(1,4)]) +
  geom_point(alpha = .5) +
  facet_wrap(~group_alpha) +
  geom_smooth(color = "black") +
  labs(x = "number coding bp within .1cM window", 
       y = "Mean mex. ancestry (by group)all in",
       main = "Effect of linked coding bp on mean mex ancestry (alpha)")
ggsave("../plots/scatter_linked_coding_bp_alpha_0.1cM_windows.png",
       height = 15, width = 20,
       units = "in", device = "png")



# plot ancestry by recombination rate
#rownames(d)<-NULL
# add in ancestry information for each population 
# at each site
d_anc <- bind_cols(d, data.frame(t(maize_anc)), 
                   data.frame(t(mex_anc))) %>%
  gather(., "pop", "alpha_mex", 17:37)
# for individual pops, plot recomb. rate vs. ancestry
d_anc %>%
  filter(., chr == "4") %>%
  filter(., pop == "maiz.Na" | pop == "maiz.Op") %>%
  ggplot(., aes(x = rate, y = alpha_mex)) +
  geom_point(cex= .1, alpha = .5, col = "black") +
  geom_smooth(col = "blue") +
  facet_wrap(~pop)
ggsave("../plots/anc_by_recomb_rate_0.1cM_windows.png", 
       height = 20, width = 25, 
       units = "in", device = "png")
# and plot coding density vs. ancestry
d_anc %>%
  #filter(., chr == "4") %>%
  #filter(., pop == "maiz.Na" | pop == "maiz.Op") %>%
  ggplot(., aes(x = coding_bp, y = alpha_mex)) +
  geom_point(cex= .1, alpha = .5, 
             col = "black") +
  geom_smooth(col = "blue") +
  facet_wrap(~pop)
ggsave("../plots/anc_by_coding_bp.png", 
       height = 20, width = 25, 
       units = "in", device = "png")
# percent coding bp
d_anc %>%
  #filter(., chr == "4") %>%
  #filter(., pop == "maiz.Na" | pop == "maiz.Op") %>%
  ggplot(., aes(x = coding_bp/width_cM, y = alpha_mex)) +
  geom_bin2d() +
  #geom_point(cex= .1, alpha = .5, 
  #           col = "black") +
  geom_smooth(col = "blue") +
  facet_wrap(~pop)
ggsave("../plots/anc_by_coding_bp_per_cM.png", 
       height = 20, width = 25, 
       units = "in", device = "png")

# plain linked bp within 0.1cM windows
p_bp <- d_anc %>%
  #filter(., chr == "4") %>%
  #filter(., pop == "maiz.Na" | pop == "maiz.Op") %>%
  ggplot(., aes(x = width_bp, y = alpha_mex)) 
#p_bp + geom_hex() # try hexbin out with smaller data
p_bp + geom_hex() +
  geom_smooth(col = "blue") +
  facet_wrap(~pop)
p_bp + geom_bin2d() +
  geom_smooth(col = "blue") +
  facet_wrap(~pop)
#p_bp + geom_point(cex= .1, alpha = .05, 
#             col = "black") +
#  geom_smooth(col = "blue") +
#  facet_wrap(~pop)
ggsave("../plots/anc_by_width_bp_0.1cM_windows_bins.png", 
       height = 20, width = 25, 
       units = "in", device = "png")


# Now I run similar analyses but with alpha and recombination and gene density statistics
# all binned in 10kb windows across the genome:
l = 10000 # length of bins
wind_10kb <- cbind(read.table("../data/refMaize/windows_10kb/gene_overlap.bed", sep = "\t", header = F, stringsAsFactors = F),
                   read.table("../data/refMaize/windows_10kb/whole_genome.recomb", header = F, stringsAsFactors = F))   
colnames(wind_10kb) <- c("chr", "start", "end", "wind_10kb_N", "n_CDS", "coding_bp", "width_bp", "perc_coding", "rate")
# which 10kb bin do variants belong in?
d_10kb <- d %>%
  #dplyr::select(chr, pos, inv_any, ZAnc_maize, ZAnc_mex, meanAlpha_maize, meanAlpha_mex) %>%
  dplyr::select(chr, pos, inv4m, meanAlpha_maize, meanAlpha_mex) %>%
  mutate(start = floor(pos/l)*l) %>%
  group_by(chr, start) %>% # any SNPs within inversion, bin is within inversion
  summarize(inv_any = Reduce("|", inv4m), 
  meanAlpha_maize = mean(meanAlpha_maize),
  meanAlpha_mex = mean(meanAlpha_mex)) %>%
  left_join(., wind_10kb, by = c("chr", "start")) %>%
  mutate(coding_density = coding_bp/(width_bp*rate/10^6*100))
d_10kb$bin_rate = with(d_10kb, cut(rate, 
                        breaks = quantile(rate, p = seq(0, 1, by = .1)),
                        right = T,
                        include.lowest = T))

# a good ~25% portion of windows this size had no estimated ancestry
# and so are dropped
1-dim(d_10kb)[1]/dim(wind_10kb)[1]
ggplot(wind_10kb, aes(x = rate, y = perc_coding)) +
  geom_point() +
  geom_smooth(color = "blue") + 
  labs(title = "across 10kb bins higher percent coding in regions with high recomb. rate",
          x = "recomb. rate - cM/Mb")
ggsave("../plots/scatter_perc_coding_by_r_10kb_windows.png",
       height = 5, width = 7,
       units = "in", device = "png")
# no qualitative difference when excluding bins with no ancestry calls
# or highlighting SNPs in inversion
d_10kb %>%
  ggplot(., aes(x = rate, y = perc_coding, color = inv_any)) +
  geom_point(alpha = .3) +
  geom_smooth(color = "blue") + 
  labs(title = "across 10kb bins WITH ANCESTRY CALLS higher percent coding in regions with high recomb. rate",
       x = "recomb. rate - cM/Mb")

# can't take true quantiles for % coding because the first several collapse into each other
# ie the windows are small enough that >50% have 0 coding bp
d_10kb$bin_perc_coding = with(d_10kb, cut(perc_coding, 
                               breaks = unique(quantile(perc_coding, p = seq(0, 1, by = .1))),
                               right = T,
                               include.lowest = T))
d_10kb$bin_coding_density = with(d_10kb, cut(coding_density, 
                                          breaks = unique(quantile(coding_density, p = seq(0, 1, by = .1))),
                                          right = T,
                                          include.lowest = T))

# plots
# violin plot alpha by recombination rate
violin_rate_alpha_10kb <- d_10kb %>%
  filter(!inv_any) %>%
  gather("group_alpha", "alpha", 
         c("meanAlpha_mex", "meanAlpha_maize")) %>%
  ggplot(., aes(y = alpha, x = bin_rate,
                      fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_violin() +
  stat_summary(fun.y=mean, geom="point", shape=3, size=2, color = "black", show.legend = F) +
  stat_summary(fun.y=median, geom="point", shape=3, size=2, color = "white", show.legend = F) +
  facet_wrap(~group_alpha) +
  labs(x = "recombination rate (10kb winds)", 
       y = "Mean mex. ancestry by group",
       title = "Effect of local recomb. on mean mex ancestry (alpha)")
plot(violin_rate_alpha_10kb)
ggsave("plots/violin_rate_alpha_10kb_windows.png",
       height = 5, width = 8,
       units = "in", device = "png")
# boxplot alpha by recombination rate
boxplot_rate_alpha_10kb <- d_10kb %>%
  filter(!inv_any) %>%
  gather("group_alpha", "alpha", 
         c("meanAlpha_mex", "meanAlpha_maize")) %>%
  ggplot(., aes(y = alpha, x = bin_rate,
                fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_boxplot() +
  facet_wrap(~group_alpha) +
  labs(x = "recombination rate (10kb winds)", 
       y = "Mean mex. ancestry by group",
       title = "Effect of local recomb. on mean mex ancestry (alpha)")
plot(boxplot_rate_alpha_10kb)
ggsave("plots/boxplot_rate_alpha_10kb_windows.png",
       height = 5, width = 8,
       units = "in", device = "png")

# plot alpha by coding density -- not true quantiles b/c so many zeros!
boxplot_perc_coding_alpha_10kb <- d_10kb %>%
  filter(!inv_any) %>%
  gather("group_alpha", "alpha", 
         c("meanAlpha_mex", "meanAlpha_maize")) %>%
  ggplot(., aes(y = alpha, x = bin_perc_coding,
                fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_boxplot() +
  facet_wrap(~group_alpha) +
  labs(x = "perc coding (10kb winds - not true quantiles)", 
       y = "Mean mex. ancestry by group",
       title = "Effect of perc coding on mean mex ancestry (alpha)")
plot(boxplot_perc_coding_alpha_10kb)
ggsave("plots/boxplot_perc_coding_alpha_10kb_windows.png",
       height = 5, width = 8,
       units = "in", device = "png")
# not super convincing but htere are so many points in low coding density side of this fig.
boxplot_coding_density_alpha_10kb <- d_10kb %>%
  filter(!inv_any) %>%
  gather("group_alpha", "alpha", 
         c("meanAlpha_mex", "meanAlpha_maize")) %>%
  ggplot(., aes(y = alpha, x = bin_coding_density,
                fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_boxplot() +
  facet_wrap(~group_alpha) +
  labs(x = "coding density bp/cM (10kb winds - not true quantiles)", 
       y = "Mean mex. ancestry by group",
       title = "Effect of coding density on mean mex ancestry (alpha)")
plot(boxplot_coding_density_alpha_10kb)
ggsave("plots/boxplot_coding_density_alpha_10kb_windows.png",
       height = 5, width = 8,
       units = "in", device = "png")

# boxplot percent coding by recombination rate
boxplot_rate_perc_coding_10kb <- d_10kb %>%
  filter(!inv_any) %>%
  ggplot(., aes(y = perc_coding, x = bin_rate)) +
  geom_boxplot() +
  labs(x = "recombination rate (10kb winds)", 
       y = "percent coding bp",
       title = "Effect of local recomb. on % coding bp")
plot(boxplot_rate_perc_coding_10kb)
ggsave("plots/boxplot_rate_perc_coding_10kb_windows.png",
       height = 5, width = 8,
       units = "in", device = "png")

# what if I bump up the scale to 100kb?
wind_100kb <- mutate(wind_10kb, 
                     start_100kb = floor(start/100000)*100000) %>%
  group_by(chr, start_100kb) %>%
  summarize(coding_bp = sum(coding_bp),
            width_bp = sum(width_bp),
            perc_coding = mean(perc_coding),
            rate = mean(rate))
d_100kb <- d %>%
  #dplyr::select(chr, pos, inv_any, ZAnc_maize, ZAnc_mex, meanAlpha_maize, meanAlpha_mex) %>%
  dplyr::select(chr, pos, inv4m, meanAlpha_maize, meanAlpha_mex) %>%
  mutate(start_100kb = floor(pos/100000)*100000) %>%
  group_by(chr, start_100kb) %>% # any SNPs within inversion, bin is within inversion
  summarize(inv_any = Reduce("|", inv4m), 
            meanAlpha_maize = mean(meanAlpha_maize),
            meanAlpha_mex = mean(meanAlpha_mex)) %>%
  left_join(., wind_100kb, 
            by = c("chr", "start_100kb")) %>%
  mutate(coding_density = coding_bp/(width_bp*rate/10^6*100))
d_100kb$bin_rate = with(d_100kb, cut(rate, 
                                   breaks = quantile(rate, p = seq(0, 1, by = .1)),
                                   right = T,
                                   include.lowest = T))
d_100kb$bin_perc_coding = with(d_100kb, cut(perc_coding, 
                                          breaks = unique(quantile(perc_coding, p = seq(0, 1, by = .1))),
                                          right = T,
                                          include.lowest = T))
d_100kb$bin_coding_density = with(d_100kb, cut(coding_density, 
                                             breaks = unique(quantile(coding_density, p = seq(0, 1, by = .1))),
                                             right = T,
                                             include.lowest = T))

d_100kb %>%
  filter(! inv_any) %>%
  ggplot(., aes(x = rate, y = perc_coding)) +
  geom_point(alpha = .3) +
  geom_smooth(color = "blue") + 
  labs(title = "across 100kb bins WITH ANCESTRY CALLS higher percent coding in regions with high recomb. rate",
       x = "recomb. rate - cM/Mb")
ggsave("plots/scatter_perc_coding_by_r_100kb_windows.png",
       height = 5, width = 7,
       units = "in", device = "png")

# boxplot alpha by recombination rate 100kb
boxplot_rate_alpha_100kb <- d_100kb %>%
  filter(!inv_any) %>%
  gather("group_alpha", "alpha", 
         c("meanAlpha_mex", "meanAlpha_maize")) %>%
  ggplot(., aes(y = alpha, x = bin_rate,
                fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_boxplot() +
  facet_wrap(~group_alpha) +
  labs(x = "recombination rate (100kb winds)", 
       y = "Mean mex. ancestry by group",
       title = "Effect of local recomb. on mean mex ancestry (alpha)")
plot(boxplot_rate_alpha_100kb)
ggsave("plots/boxplot_rate_alpha_100kb_windows.png",
       height = 5, width = 8,
       units = "in", device = "png")

# plot alpha by % coding -- not true quantiles b/c so many zeros!
boxplot_perc_coding_alpha_100kb <- d_100kb %>%
  filter(!inv_any) %>%
  gather("group_alpha", "alpha", 
         c("meanAlpha_mex", "meanAlpha_maize")) %>%
  ggplot(., aes(y = alpha, x = bin_perc_coding,
                fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_boxplot() +
  facet_wrap(~group_alpha) +
  labs(x = "perc coding (100kb winds - not true quantiles)", 
       y = "Mean mex. ancestry by group",
       title = "Effect of perc coding on mean mex ancestry (alpha)")
plot(boxplot_perc_coding_alpha_100kb)
ggsave("plots/boxplot_perc_coding_alpha_100kb_windows.png",
       height = 5, width = 8,
       units = "in", device = "png")

# plot alpha by coding density -- not true quantiles b/c so many zeros!
boxplot_coding_density_alpha_100kb <- d_100kb %>%
  filter(!inv_any) %>%
  gather("group_alpha", "alpha", 
         c("meanAlpha_mex", "meanAlpha_maize")) %>%
  ggplot(., aes(y = alpha, x = bin_coding_density,
                fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_boxplot() +
  facet_wrap(~group_alpha) +
  labs(x = "coding density bp/cM (100kb winds - not true quantiles)", 
       y = "Mean mex. ancestry by group",
       title = "Effect of coding density on mean mex ancestry (alpha)")
plot(boxplot_coding_density_alpha_100kb)
ggsave("plots/boxplot_coding_density_alpha_100kb_windows.png",
       height = 5, width = 8,
       units = "in", device = "png")


# what about 1Mb?
wind_1000kb <- mutate(wind_10kb, 
                     start_1000kb = floor(start/1000000)*1000000) %>%
  group_by(chr, start_1000kb) %>%
  summarize(coding_bp = sum(coding_bp),
            width_bp = sum(width_bp),
            perc_coding = mean(perc_coding),
            rate = mean(rate))
d_1000kb <- d %>%
  dplyr::select(chr, pos, inv_any, ZAnc_maize, ZAnc_mex, meanAlpha_maize, meanAlpha_mex) %>%
  mutate(start_1000kb = floor(pos/1000000)*1000000) %>%
  group_by(chr, start_1000kb) %>% # any SNPs within inversion, bin is within inversion
  summarize(inv_any = Reduce("|", inv_any), 
            meanAlpha_maize = mean(meanAlpha_maize),
            meanAlpha_mex = mean(meanAlpha_mex)) %>%
  left_join(., wind_1000kb, 
            by = c("chr", "start_1000kb")) %>%
  mutate(coding_density = coding_bp/(width_bp*rate/10^6*100))
d_1000kb$bin_rate = with(d_1000kb, cut(rate, 
                                     breaks = quantile(rate, p = seq(0, 1, by = .1)),
                                     right = T,
                                     include.lowest = T))
d_1000kb$bin_perc_coding = with(d_1000kb, cut(perc_coding, 
                                            breaks = unique(quantile(perc_coding, p = seq(0, 1, by = .1))),
                                            right = T,
                                            include.lowest = T))
d_1000kb$bin_coding_density = with(d_1000kb, cut(coding_density, 
                                               breaks = unique(quantile(coding_density, p = seq(0, 1, by = .1))),
                                               right = T,
                                               include.lowest = T))

d_1000kb %>%
  filter(! inv_any) %>%
  ggplot(., aes(x = rate, y = perc_coding)) +
  geom_point(alpha = .3) +
  geom_smooth(color = "blue") + 
  labs(title = "across 1000kb bins WITH ANCESTRY CALLS higher percent coding in regions with high recomb. rate",
       x = "recomb. rate - cM/Mb")
ggsave("../plots/scatter_perc_coding_by_r_1000kb_windows.png",
       height = 5, width = 7,
       units = "in", device = "png")

# boxplot alpha by recombination rate 1000kb
boxplot_rate_alpha_1000kb <- d_1000kb %>%
  filter(!inv_any) %>%
  gather("group_alpha", "alpha", 
         c("meanAlpha_mex", "meanAlpha_maize")) %>%
  ggplot(., aes(y = alpha, x = bin_rate,
                fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_boxplot() +
  facet_wrap(~group_alpha) +
  labs(x = "recombination rate (1000kb winds)", 
       y = "Mean mex. ancestry by group",
       title = "Effect of local recomb. on mean mex ancestry (alpha)")
plot(boxplot_rate_alpha_1000kb)
ggsave("../plots/boxplot_rate_alpha_1000kb_windows.png",
       height = 5, width = 8,
       units = "in", device = "png")

# plot alpha by % coding
boxplot_perc_coding_alpha_1000kb <- d_1000kb %>%
  filter(!inv_any) %>%
  gather("group_alpha", "alpha", 
         c("meanAlpha_mex", "meanAlpha_maize")) %>%
  ggplot(., aes(y = alpha, x = bin_perc_coding,
                fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_boxplot() +
  facet_wrap(~group_alpha) +
  labs(x = "perc coding (1000kb winds)", 
       y = "Mean mex. ancestry by group",
       title = "Effect of perc coding on mean mex ancestry (alpha)")
plot(boxplot_perc_coding_alpha_1000kb)
ggsave("../plots/boxplot_perc_coding_alpha_1000kb_windows.png",
       height = 5, width = 8,
       units = "in", device = "png")

# plot alpha by coding density -- not true quantiles b/c so many zeros!
boxplot_coding_density_alpha_1000kb <- d_1000kb %>%
  filter(!inv_any) %>%
  gather("group_alpha", "alpha", 
         c("meanAlpha_mex", "meanAlpha_maize")) %>%
  ggplot(., aes(y = alpha, x = bin_coding_density,
                fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_boxplot() +
  facet_wrap(~group_alpha) +
  labs(x = "coding density bp/cM (1000kb winds)", 
       y = "Mean mex. ancestry by group",
       title = "Effect of coding density on mean mex ancestry (alpha)")
plot(boxplot_coding_density_alpha_1000kb)
ggsave("../plots/boxplot_coding_density_alpha_1000kb_windows.png",
       height = 5, width = 8,
       units = "in", device = "png")


# plot mean population ancestry across the genome:
maize_anc <- all_anc[ , maize_pops]
mexicana_anc <- all_anc[ , mexicana_pops]
d %>%
  ggplot(aes(x = pos, y = meanAlpha_maize, color = inv4m)) +
  geom_point(size = .1) +
  facet_wrap(~chr) + 
  ylab("mean mexicana freq. across all individuals") +  
  ggtitle("frequency of mexicana ancestry in maize") +
  geom_abline(intercept = mean(d$meanAlpha_maize), slope = 0, color = "darkgrey")
ggsave("plots/mean_mexicana_in_maize_individuals_genomewide.png",
       height = 10, width = 14,
       units = "in", device = "png")
d %>%
  ggplot(aes(x = pos, y = pop_meanAlpha_maize, color = inv4m)) +
  geom_point(size = .1) +
  facet_wrap(~chr) + 
  ylab("mean mexicana freq. across all populations") +
  ggtitle("frequency of mexicana ancestry in maize") +
  geom_abline(intercept = mean(d$pop_meanAlpha_maize), slope = 0, color = "darkgrey")
ggsave("plots/mean_mexicana_in_maize_pops_genomewide.png",
       height = 10, width = 14,
       units = "in", device = "png")
summary(d$meanAlpha_maize) # mean is around 31% mexicana frequency

# the two plots above don't look very different because it's taking a mean of populations or pooling
# across populations makes little difference -- the populations are sampled roughly evenly & track together

# again little difference if I look at pooled individuals or pooled pops to take a mean ancestry
d %>%
  ggplot(aes(x = meanAlpha_maize, y = pop_meanAlpha_maize, color = chr)) +
  geom_point(size = .1) + 
  ggtitle("population vs. indivudal mean mexicana ancestry in maize")
ggsave("plots/pop_vs_ind_mean_mexicana_ancestry_in_maize.png",
       height = 6, width = 6,
       units = "in", device = "png")


# what does the pattern of mexicana ancestry look like in sympatric mexicana?
d %>%
  ggplot(aes(x = meanAlpha_mex, y = pop_meanAlpha_mex, color = chr)) +
  geom_point(size = .1) + 
  ggtitle("population vs. indivudal mean mexicana ancestry in maize")
ggsave("plots/pop_vs_ind_mean_mexicana_ancestry_in_mexicana.png",
       height = 6, width = 6,
       units = "in", device = "png")
# mexicana ancestry across the genome in mexicana pops
d %>%
  ggplot(aes(x = pos, y = meanAlpha_mex, color = inv4m)) +
  geom_point(size = .1) +
  facet_wrap(~chr) + 
  ylab("mean mexicana freq. across all individuals") +  
  ggtitle("frequency of mexicana ancestry in mexicana") +
  geom_abline(intercept = mean(d$meanAlpha_mex), slope = 0, color = "darkgrey")
ggsave("plots/mean_mexicana_in_mexicana_individuals_genomewide.png",
       height = 10, width = 14,
       units = "in", device = "png")

# what is the correlation between mexicana ancestry in maize vs. in mexicana?
d %>%
  arrange(., inv4m) %>% # put on top to make more visible
  ggplot(aes(x = meanAlpha_maize, y = meanAlpha_mex, color = inv4m)) +
  geom_point(size = .01) + 
  ggtitle(paste0("mean mexicana ancestry in maize vs. in mexicana, cor=", cor(d$meanAlpha_maize, d$meanAlpha_mex)))
# some of this is positively selected alleles from 1 ancestry, e.g. inv4m
# I can check if this correlation goes down when I restrict to higher coverage individuals 
# (ie. is it partly caused by inability to detect a particular ancestry in some parts of the genome)
ggsave("plots/mean_mexicana_ancestry_in_maize_vs_in_mexicana.png",
       height = 6, width = 6,
       units = "in", device = "png")

# what does the picture look like when I split it out by population?
d.pops <- bind_cols(d, all_anc) %>%
  gather(., "pop", "meanAnc", c("pop_meanAlpha_maize", "pop_meanAlpha_mex", "meanAlpha_maize", "meanAlpha_mex", colnames(all_anc)))
d.pops %>%
  filter(pop %in% maize_pops) %>%
  ggplot(aes(x = pos, y = meanAnc, color = pop)) +
  geom_point(size = .1) +
  facet_wrap(~chr) +
  ggtitle("mean mexicana ancestry in diff. maize populations")
ggsave("plots/mexicana_ancestry_in_maize_pops_genomewide.png",
       height = 10, width = 14,
       units = "in", device = "png")

# can I see outliers by elevation?
d.pops %>%
  filter(pop %in% maize_pops) %>%
  left_join(., meta.pops, by = "pop") %>%
  ggplot(aes(x = pos, y = meanAnc, color = ELEVATION)) +
  geom_point(size = .1, alpha = .2) +
  facet_wrap(~chr) +
  ylab("mean population mexicana ancestry") +
  ggtitle("mean mexicana ancestry in maize populations by elevation") +
  scale_color_gradientn(colors = brewer.pal(7, "YlGnBu")) # change the colors
ggsave("plots/mexicana_ancestry_in_maize_pops_by_elevation_genomewide.png",
       height = 10, width = 14,
       units = "in", device = "png")


d.pops %>%
  filter(pop %in% maize_pops) %>%
  filter(chr == 4) %>% # just chr4
  left_join(., meta.pops, by = "pop") %>%
  ggplot(aes(x = pos, y = meanAnc, color = ELEVATION)) +
  geom_point(size = .1, alpha = .2) +
  facet_wrap(~chr) +
  ylab("mean population mexicana ancestry") +
  ggtitle("mean mexicana ancestry in maize populations by elevation") +
  scale_color_gradientn(colors = brewer.pal(7, "YlGnBu")) # change the colors
ggsave("plots/mexicana_ancestry_in_maize_pops_by_elevation_chr4.png",
       height = 10, width = 14,
       units = "in", device = "png")

d.pops %>%
  filter(pop %in% maize_pops | pop %in% mexicana_pops) %>%
  left_join(., meta.pops, by = "pop") %>%
  group_by(inv4m, pop, ELEVATION, LOCALITY, zea) %>%
  summarise(mex_anc_freq = mean(meanAnc)) %>%
  ggplot(aes(x = ELEVATION, y = mex_anc_freq, color = inv4m)) +
  geom_smooth(method = "lm") +
  geom_point() +
  facet_wrap(~zea) +
  ylab("mean population mexicana ancestry") +
  ggtitle("mean mexicana ancestry in and out of inversion")
ggsave("plots/lm_mexicana_ancestry_inv4m_vs_genomewide_background_by_elevation.png",
       height = 10, width = 14,
       units = "in", device = "png")

# can I get a PCA of just the inversion?
# plot ZAnc ~ environment across the genome?
# calculate ZBeta. Does anything have negative ZBeta? 
# (neg. association between mexicana anc and elevation should be rare compared to pos ZBeta)

# helper file with some useful functions
source("../../covAncestry/forqs_sim/k_matrix.R") # import useful functions
# function for calculating statistic:
make_calcs = function(pop_anc){
  # calculate mean population frequency across all snps
  pop_alpha = apply(pop_anc, 1, mean)
  # calculate K matrix for populations
  pop_K = calcK(ancFreqMatrix = pop_anc, alpha = pop_alpha)
  pop_InvL = calcInvL(pop_K)
  # for each locus, calculate ZAnc statistic
  pop_ZAnc = apply(pop_anc, 2, function(l) ZAnc(ancFreq = l, invL = pop_InvL, alpha = pop_alpha))
  pop_p_ZAnc = calcProbZAnc(ZAnc = pop_ZAnc, nPop = length(pop_alpha), logProb = F)
  return(list(alpha = pop_alpha, K = pop_K, InvL = pop_InvL, ZAnc = pop_ZAnc, p_ZAnc = pop_p_ZAnc))
}
zAnc_mexicana = make_calcs(t(mexicana_anc))
zAnc_maize = make_calcs(t(maize_anc))
meta.pops$pop[meta.pops$zea == "maize"] == colnames(maize_anc)
zBeta_elev_maize <- apply(maize_anc, 1, function(l) # calculate slope for each locus
                          #zBeta(ancFreq = t(maize_anc[1,]),
  zBeta(ancFreq = l,
        envWeights = meta.pops$ELEVATION[meta.pops$zea == "maize"], 
        invL = zAnc_maize$InvL, 
        alpha = zAnc_maize$alpha))
d$zBeta_elev_maize <- zBeta_elev_maize
hist(d$zBeta_elev_maize)

# returning more information about the model
zBeta2_elev_maize <- apply(maize_anc, 1, function(l) # calculate slope for each locus
  #zBeta(ancFreq = t(maize_anc[1,]),
  zBeta2(ancFreq = l,
        envWeights = meta.pops$ELEVATION[meta.pops$zea == "maize"], 
        invL = zAnc_maize$InvL, 
        alpha = zAnc_maize$alpha))
zb2 <- data.frame(t(zBeta2_elev_maize)) # data frame is easier to work with
with(zb2, plot(r_squared ~ slope.zEnv, col = ifelse(d$inv4m, "blue", "black")))
with(zb2, plot(sum_sq_res ~ slope.zEnv, col = ifelse(d$inv4m, "blue", "black")))

# not centering environment:
zBeta3_elev_maize <- apply(maize_anc, 1, function(l) # calculate slope for each locus
  zBeta3(ancFreq = l,
         envWeights = meta.pops$ELEVATION[meta.pops$zea == "maize"], 
         invL = zAnc_maize$InvL, 
         alpha = zAnc_maize$alpha)) 
zb3_elev <- data.frame(t(zBeta3_elev_maize)) # data frame is easier to work with
zBeta3_highmex_maize <- apply(maize_anc, 1, function(l) # calculate slope for each locus
  zBeta3(ancFreq = l,
         envWeights = rep(1, 14), 
         invL = zAnc_maize$InvL, 
         alpha = zAnc_maize$alpha)) 
zb3_hmex <- data.frame(t(zBeta3_highmex_maize)) # data frame is easier to work with

with(zb3_elev, plot(r_squared ~ zEnv, col = ifelse(d$inv4m, "blue", "black"), main = "zBeta3 association mexicana ancestry with high elevation in maize"))
with(zb3_elev, plot(sum_sq_res ~ zEnv, col = ifelse(d$inv4m, "blue", "black"), main = "zBeta3 association mexicana ancestry with high elevation in maize"))
with(zb3_hmex, plot(r_squared ~ zEnv, col = ifelse(d$inv4m, "blue", "black"), main = "zBeta3 universally high mexicana anc in maize"))
with(zb3_hmex, plot(sum_sq_res ~ zEnv, col = ifelse(d$inv4m, "blue", "black"), main = "zBeta3 universally high mexicana anc in maize"))

zBeta3_elev_int_maize <- apply(maize_anc, 1, function(l)
  zBeta3_int(ancFreq = l,
             envWeights = meta.pops$ELEVATION[meta.pops$zea == "maize"],
             invL = zAnc_maize$InvL,
             alpha = zAnc_maize$alpha))
zb3_elev_int <- data.frame(t(zBeta3_elev_int_maize))


# plot outliers against their population means
d %>%
  bind_cols(., zb3_hmex) %>%
  filter(., pval < .01) %>%
  ggplot(aes(x = pop_meanAlpha_maize, y = zEnv, color = log(pval))) +
  geom_point(size = .1) +
  ggtitle("zAnc outliers, p < .01")
ggsave(paste0("plots/zAnc_top_mex_outliers_in_maize_1percent.png"),
       device = "png",
       width = 10, height = 8, units = "in")
d %>%
  bind_cols(., zb3_hmex) %>%
  ggplot(aes(x = pop_meanAlpha_maize, y = zEnv, color = log(pval))) +
  geom_point(size = .1) +
  ggtitle("zAnc vs. population mean mexicana")
ggsave(paste0("plots/zAnc_pop_mean_mexicana_in_maize_allLoci.png"),
       device = "png",
       width = 10, height = 8, units = "in")

d %>% # color by sum of squared residuals. okay actually that's only relevant as a relative value compared to the null model residuals
  bind_cols(., zb3_hmex) %>%
  filter(., r_squared > .5) %>%
  ggplot(aes(x = pop_meanAlpha_maize, y = zEnv, color = r_squared)) +
  geom_point(size = .1)

# calculate zAnc at each locus. 
zAnc3 <- t(apply(maize_anc, 1, function(l) zAnc(ancFreq = l, 
                                      invL = zAnc_maize$InvL, 
                                      alpha = zAnc_maize$alpha)))

# what do individual top outliers look like?
# e.g. frequency of pval < .05 zElev outliers across pops, colored by inversion

bind_cols(maize_anc, d) %>%
  filter(., zb3_elev$pval < .05) %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = mex_freq, fill = inv4m)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# now what if I let variation with environment have an intercept?
bind_cols(maize_anc, d) %>%
  bind_cols(., zb3_elev_int) %>%
  filter(., pval_zEnv < .01) %>%
  .[c(T, rep(F, 100)), ] %>%
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = mex_freq, shape = inv4m, color = zEnv, group = pos)) +
  geom_point() +
  geom_line() +
  facet_wrap(~chr) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset low p-val outliers from zElev test - with rotated intercept zInt")
ggsave(paste0("plots/a_few_example_loci_zElev_outliers_rotated_intercept.png"),
       device = "png",
       width = 10, height = 8, units = "in")

bind_cols(maize_anc, d) %>%
  bind_cols(., zb3_elev) %>%
  filter(., pval < .01) %>%
  .[ c(T, rep(F, 240)), ] %>% # take subset
  gather(., "pop", "mex_freq", maize_pops) %>%
  left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = mex_freq, shape = inv4m, color = zEnv, group = pos)) +
  geom_point() + 
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("subset low p-val outliers from zElev test - no intercept")
ggsave(paste0("plots/a_few_example_loci_zElev_outliers_no_intercept.png"),
       device = "png",
       width = 10, height = 8, units = "in")


# what are the dominant patterns in my data?
# note: default algorithm is Hartigan-Wong. My other choices include 'Lloyd' and 'MacQueen'
set.seed(10) # set seed for k-means clustering
kmeans(zAnc3, 1)$centers # good, with k=1 one cluster, the mean should be all ~ zeros
k1 <- kmeans(zAnc3, centers = 1, nstart = 25, iter.max = 25)
k2 <- kmeans(zAnc3, centers = 2, nstart = 25, iter.max = 25)
k2$centers
k3 <- kmeans(zAnc3, centers = 3, nstart = 25, iter.max = 50)
#k4 <- kmeans(zAnc3, centers = 4, nstart = 25, iter.max = 50, trace = T) # produces a lot of output
k4 <- kmeans(zAnc3, centers = 4, nstart = 25, iter.max = 50)
k5 <- kmeans(zAnc3, centers = 5, nstart = 25, iter.max = 50)
k10 <- kmeans(zAnc3, centers = 10, nstart = 25, iter.max = 50) # 22 warnings, but 25 starts? it's choosing the best one out of the ones that don't fail I assume
# I increased iterations and number of starts to avoid nonconvergence errors:
# 16: Quick-TRANSfer stage steps exceeded maximum (= 19105350)

kmeans(maize_anc, 1)$centers

# first look at clusters in all positions (not just outliers):
colnames(maize_anc) <- maize_pops
maize_anc %>%
  mutate(k.2 = as.factor(k2$cluster)) %>% # which cluster does this locus belong to?
  gather(., "pop", "mex_freq", maize_pops) %>%
  ggplot(aes(x = pop, y = mex_freq, fill = k.2)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
maize_anc %>% # too big a plot to view properly
  mutate(k.1 = as.factor(k1$cluster)) %>%
  mutate(k.2 = as.factor(k2$cluster)) %>%
  mutate(k.3 = as.factor(k3$cluster)) %>%
  mutate(k.4 = as.factor(k4$cluster)) %>%
  mutate(k.5 = as.factor(k5$cluster)) %>%
  mutate(k.10 = as.factor(k10$cluster)) %>% 
  gather(., "pop", "mex_freq", maize_pops) %>%
  gather(., "k", "cluster", paste("k", c(1:5, 10), sep = ".")) %>%
  ggplot(aes(x = pop, y = mex_freq, fill = cluster)) +
  geom_boxplot() +
  facet_wrap(~k) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("all loci maize - k-means clusters")
# save individual plots
ks <- list(k1, k2, k3, k4, k5, k10)
ks_names <- paste0("k", c(1:5, 10))
plot_ks <- lapply(1:length(ks), function(i)
  maize_anc %>%
    mutate(cluster = as.factor(ks[[i]]$cluster)) %>% # which cluster does this locus belong to?
    gather(., "pop", "mex_freq", maize_pops) %>%
    left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
    ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = mex_freq, fill = cluster)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(paste("k-means clustering:", ks_names[i])))
for (i in 1:length(ks)){
  plot(plot_ks[[i]])
  ggsave(paste0("plots/k_means_all_loci_maize_", ks_names[i], ".png"),
         device = "png",
         width = 10, height = 8, units = "in")
}
  
 
# Now get zTz -- basically the sum of residual errors to the null model
# zTz
zTz3 <- apply(zAnc3, 1, function(i) t(i) %*% i)
d %>%
  mutate(zTz = zTz3) %>%
  bind_cols(., zb3_hmex) %>%
  filter(zTz >= quantile(zTz, .95)) %>% # limit to top 1%
  mutate(log_pval = log(pval)) %>%
  ggplot(aes(x = meanAlpha_maize, y = zEnv, color = log_pval)) +
  geom_point(size = .1) +
  ggtitle("only showing top 5% zTz outliers; universal selection model")
d %>% # note: inv4m is not in top 1% of zTz but it is in top 5% of zTz
  mutate(zTz = zTz3) %>%
  bind_cols(., zb3_elev) %>%
  filter(zTz >= quantile(zTz, .95)) %>% # limit to top 1%
  mutate(log_pval = log(pval)) %>%
  ggplot(aes(x = meanAlpha_maize, y = zEnv, color = log_pval)) +
  #ggplot(aes(x = meanAlpha_maize, y = zEnv, color = inv4m)) +
  geom_point(size = .1) +
  ggtitle("only showing top 5% zTz outliers; elevation gradient model")

# now cluster just within the top 1% outliers
top1.zAnc3 <- zAnc3[zTz3 >= quantile(zTz3, .99), ]
top1.ks <- lapply(1:5, function(k)
  kmeans(top1.zAnc3, centers = k, nstart = 25, iter.max = 50))
top1.ks_names <- paste0("k", 1:5)
top1.plot_ks <- lapply(1:length(top1.ks), function(i)
  maize_anc[zTz3 >= quantile(zTz3, .99), ] %>%
    mutate(cluster = as.factor(top1.ks[[i]]$cluster)) %>% # which cluster does this locus belong to?
    gather(., "pop", "mex_freq", maize_pops) %>%
    left_join(., meta.pops[ , c("pop", "ELEVATION", "LOCALITY")], by = "pop") %>%
    ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = mex_freq, fill = cluster)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(paste("top 1% zTz hits: k-means clustering:", top1.ks_names[i])))
for (i in 1:length(top1.ks)){
  plot(top1.plot_ks[[i]])
  ggsave(paste0("plots/k_means_top1percent_loci_maize_", top1.ks_names[i], ".png"),
         device = "png",
         width = 10, height = 8, units = "in")
}



# what does zEnv look like genomewide for elevation?
d %>%
  mutate(zTz = zTz3) %>%
  bind_cols(., zb3_hmex) %>%
  ggplot(aes(x = pos, y = zEnv, color = zTz, shape = inv4m)) +
  geom_point(size = .1) +
  ggtitle("zElev outliers across the genome") +
  facet_wrap(~chr)
d %>%
  mutate(zTz = zTz3) %>%
  bind_cols(., zb3_hmex) %>%
  mutate(., log_pval = log(pval)) %>%
  ggplot(aes(x = pos, y = zEnv, color = zTz)) +
  geom_point(size = .1) +
  ggtitle("zAnc outliers across the genome") +
  facet_wrap(~chr)
ggsave("plots/zAnc_in_maize_pops_genomewide.png",
       height = 10, width = 14,
       units = "in", device = "png")
d %>%
  mutate(zTz = zTz3) %>%
  bind_cols(., zb3_elev) %>%
  mutate(., log_pval = log(pval)) %>%
  ggplot(aes(x = pos, y = zEnv, color = zTz)) +
  geom_point(size = .1) +
  ggtitle("zElev outliers across the genome") +
  facet_wrap(~chr)
ggsave("plots/zElev_in_maize_pops_genomewide.png",
       height = 10, width = 14,
       units = "in", device = "png")

# with intercept zInt for elevation:
d %>%
  mutate(zTz = zTz3) %>%
  bind_cols(., zb3_elev_int) %>%
  mutate(., log_pval = log(pval_zEnv)) %>%
  ggplot(aes(x = pos, y = zEnv, color = (inv4m | inv9e_or_d | inv9d_or_e))) +
  geom_point(size = .1) + 
  ggtitle("zElev outliers across the genome (w/ intercept)") +
  facet_wrap(~chr) # color inversion on chr9 too .. do I have the correct coordinates for these? is the length reasonable? are there others?
ggsave("plots/zElev_in_maize_pops_genomewide_w_intercept.png",
       height = 10, width = 14, 
       units = "in", device = "png")

# what is the qualitative diff of centering? zb2 vs. zb3:
plot(zb2$slope.zEnv ~ zb3_elev$slope.zEnv) # eh blob
cor(zb2$slope.zEnv, zb3_elev$slope.zEnv) # very weak correlation, or = .04
plot(zb2$pval_slope ~ zb3_elev$pval_slope) # black box (everywhere)
cor(zb2$pval_slope, zb3_elev$pval_slope) # slightly negative, cor -.06
plot(zb2$sum_sq_res ~ zb3_elev$sum_sq_res) # highly correlated
cor(zb2$sum_sq_res, zb3_elev$sum_sq_res) # cor = .96
cor(zb3_elev_int$zEnv, zb3_elev$zEnv) # poorly correlated
# plot zb3_elev_int$zEnv vs. simple beta elevation


# test set: there are 14 populations of maize. I will create some test sets with 14 pops
# that have perfect environmental correlations
# each row is a made-up ancestry frequency observation; each column is a population 
test_case_names <- c("all_1", "last_pop_only", "sixth_pop_only", "first_pop_only", "linear_env", "first_one_out", "last_one_out", "linear_env+100", "linear_env*2", "linear_env*2+100")

test_anc <- rbind(rep(1, 14), # all selected for mexicana ancestry
                  c(zAnc_maize$alpha[1:13], 1), # only last population selected
                  c(zAnc_maize$alpha[1:5], 1, zAnc_maize$alpha[7:14]), # 6th population selected
                  c(1, zAnc_maize$alpha[2:14]), # last population selected
                  seq(0, 1, length.out = 14),
                  c(0, rep(1, 13)),
                  c(rep(1,13), 0),
                  seq(0, 1, length.out = 14) + 100,
                  seq(0, 1, length.out = 14)*2,
                  seq(0, 1, length.out = 14)*2 + 100) # linear association with simple linear environmental gradient: e.g. seq(100, 200, length.out = 14)
                  # could also test false positive cases: e.g. one outlier but no real other association with

# I'll use my test_anc cases with the empirical K matrix:
test_env <- rbind(rep(1, 14), # all selected for mexicana ancestry
                   c(zAnc_maize$alpha[1:13], 1), # only last population selected
                   c(zAnc_maize$alpha[1:5], 1, zAnc_maize$alpha[7:14]), # 6th population selected
                   c(1, zAnc_maize$alpha[2:14]), # last population selected
                   seq(0, 1, length.out = 14))

# transformed test environments:
apply(test_env, 1, function(e) zAnc_maize$InvL %*% e)
apply(test_env, 1, function(e) zAnc_maize$InvL %*% (e - mean(e)))

zTest2 <- lapply(1:nrow(test_env), function(e)
  apply(test_anc, 1, function(l) # calculate slope for each locus
  #zBeta(ancFreq = t(maize_anc[1,]),
  zBeta2(ancFreq = l,
         envWeights = test_env[e, ], 
         invL = zAnc_maize$InvL, 
         alpha = zAnc_maize$alpha)))

zTest3 <- lapply(1:nrow(test_env), function(e)
  apply(test_anc, 1, function(l) # calculate slope for each locus
    #zBeta(ancFreq = t(maize_anc[1,]),
    zBeta3(ancFreq = l,
           envWeights = test_env[e, ], 
           invL = zAnc_maize$InvL, 
           alpha = zAnc_maize$alpha)))
a2 <- do.call(rbind,
             lapply(1:length(zTest2), function(l)
  data.frame(t(zTest2[[l]])) %>%
  mutate(test_model = test_case_names) %>%
  mutate(true_model = test_case_names[l])))
a2 %>%
  ggplot(aes(x = test_model, y = r_squared, color = true_model, shape = (true_model == test_model))) +
  geom_point() + 
  ggtitle("mean centered environment test cases")
a2 %>%
  gather(., "model_output", "value", c("intercept..Intercept.", "slope.zEnv", "r_squared")) %>%
  ggplot(aes(x = test_model, y = value, color = true_model, shape = (true_model == test_model))) +
  geom_point() + 
  ggtitle("mean centered environment test cases") +
  facet_wrap(~model_output) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/test_cases_environment_ZAnc_model_values_centered_env.png",
       height = 10, width = 10, device = "png", units = "in")

a2 %>%
  gather(., "model_output", "value", c("sum_sq_res", "pval_intercept", "pval_slope.zEnv")) %>%
  mutate(., log_value = log(value)) %>%
  ggplot(aes(x = test_model, y = log_value, color = true_model, shape = (true_model == test_model))) +
  geom_point() + 
  ggtitle("mean centered environment test cases") +
  facet_wrap(~model_output) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/test_cases_environment_ZAnc_fit_values_centered_env.png",
       height = 10, width = 10, device = "png", units = "in")


a3 <- do.call(rbind,
              lapply(1:length(zTest3), function(l)
                data.frame(t(zTest3[[l]])) %>%
                  mutate(test_model = test_case_names) %>%
                  mutate(true_model = test_case_names[l])))
a3 %>%
  ggplot(aes(x = test_model, y = r_squared, color = true_model, shape = (true_model == test_model))) +
  geom_point() + 
  ggtitle("not mean centered environment test cases")
a3 %>%
  gather(., "model_output", "value", c("zEnv", "r_squared")) %>%
  mutate(., log_value = log(value)) %>%
  ggplot(aes(x = true_model, y = log_value, color = test_model, shape = (true_model == test_model))) +
  geom_point() + 
  ggtitle("not mean centered environment test cases") +
  facet_wrap(~model_output) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/test_cases_environment_ZAnc_model_values.png",
       height = 10, width = 10, device = "png", units = "in")

a3 %>%
  gather(., "model_output", "value", c("sum_sq_res", "pval")) %>%
  mutate(., log_value = log(value)) %>%
  ggplot(aes(x = true_model, y = log_value, color = test_model, shape = (true_model == test_model))) +
  geom_point() + 
  ggtitle("not mean centered environment test cases") +
  facet_wrap(~model_output) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/test_cases_environment_ZAnc_fit_values.png",
       height = 10, width = 10, device = "png", units = "in")

# add in an intercept of 1's vectors for the environmental correlation case
# b/c I don't want the environmental effect to be constrained to the non-centered environment
# but I won't test the c(1,1,1,1) model with this
zTest3_int <- lapply(1:nrow(test_env), function(e)
  apply(test_anc[ , ], 1, function(l) # calculate slope for each locus
    #zBeta(ancFreq = t(maize_anc[1,]),
    zBeta3_int(ancFreq = l,
           envWeights = test_env[e, ], 
           invL = zAnc_maize$InvL, 
           alpha = zAnc_maize$alpha)))


# get full linear models as output
# testing just the linear environment true model
zTest3_int_lm <- lapply(1:nrow(test_env), function(e)
  apply(test_anc, 1, function(l) # calculate slope for each locus
    #zBeta(ancFreq = t(maize_anc[1,]),
    zBeta3_int_lm(ancFreq = l,
               envWeights = test_env[e, ], 
               invL = zAnc_maize$InvL, 
               alpha = zAnc_maize$alpha)))
lapply(zTest3_int_lm, function(m) summary(m))
length(zTest3_int_lm)


# hmm..problem is that r-squared doesn't have the same meaning if there is an intercept or not
# (all models are being compared to free intercept model, not zTz now)
a3_int <- do.call(rbind,
              lapply(1:length(zTest3_int), function(l)
                data.frame(t(zTest3_int[[l]])) %>%
                  mutate(test_model = test_case_names) %>%
                  mutate(true_model = test_case_names[l])))
a3_int %>%
  gather(., "model_output", "value", c("zEnv", "zInt", "r_squared")) %>%
  #mutate(., log_value = log(value)) %>%
  ggplot(aes(x = true_model, y = value, color = test_model, shape = (true_model == test_model))) +
  geom_point() + 
  ggtitle("not mean centered environment test cases + intercept") +
  facet_wrap(~model_output, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/test_cases_environment_ZAnc_model_values_3int.png",
       height = 10, width = 10, device = "png", units = "in")

a3_int %>%
  gather(., "model_output", "value", c("sum_sq_res", "pval_zEnv")) %>%
  mutate(., log_value = log(value)) %>%
  ggplot(aes(x = true_model, y = log_value, color = test_model, shape = (true_model == test_model))) +
  geom_point() + 
  ggtitle("not mean centered environment test cases + intercept") +
  facet_wrap(~model_output) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/test_cases_environment_ZAnc_fit_values_3int.png",
       height = 10, width = 10, device = "png", units = "in")



# instead of calculating zTz in it's own framework -- can I just put in a vector of 0's for environment?


# straight regression, no rotation:
ancBeta_elev_maize <- apply(maize_anc, 1, function(l) # calculate slope for each locus
  ancBeta(ancFreq = l,
        envWeights = meta.pops$ELEVATION[meta.pops$zea == "maize"], 
        alpha = zAnc_maize$alpha))
d$ancBeta_elev_maize <- ancBeta_elev_maize
hist(d$ancBeta_elev_maize)

# simple regression, no transformation at all:
simpleBeta = function(ancFreq, envWeights){
  # calculates the association between untransformed ancestry frequencies
  # and untransformed environmental weights (e.g. c(1,0,1,0) if discrete) or c(altitude1, altitude2) if continuous
  # where the environment is a selective environment (or apriori assumed relative strength of selection within the set of environments)
  if (!(length(envWeights == length(ancFreq)))){
    stop("Oops! dimensions for envWeights, ancFreq and alpha must match")
  }
  devEnv = (envWeights - mean(envWeights)) # mean center environment
  fit = lm(ancFreq ~ devEnv)
  # summary(fit)
  fit$coefficients["devEnv"] # output slope of simple regression
}

simpleBeta_elev_maize <- apply(maize_anc, 1, function(l) # calculate slope for each locus
  simpleBeta(ancFreq = l,
          envWeights = meta.pops$ELEVATION[meta.pops$zea == "maize"]))
d$simpleBeta_elev_maize <- simpleBeta_elev_maize

# finding outliers by squaring difference from expected (expected = 0)
d$ztz <- apply(maize_anc, 1, function(l) # calculate locus by locus
  sum(zAnc(ancFreq = l,
          invL = zAnc_maize$InvL, 
          alpha = zAnc_maize$alpha)^2))

d %>%
  ggplot(aes(x = pos, y = zBeta_elev_maize, color = inv4m)) +
  geom_point(size = .1, alpha = .2) +
  facet_wrap(~chr) +
  ylab("slope zAnc mexicana ancestry ~ elevation") +
  ggtitle("Where in the maize genome is mexicana ancestry best predicted by high elevation?")
ggsave("plots/mexicana_ancestry_in_maize_pops_slope_regression_elevation_wholegenome.png",
       height = 10, width = 14,
       units = "in", device = "png")

d %>%
  ggplot(aes(x = pos, y = ancBeta_elev_maize, color = inv4m)) +
  geom_point(size = .1, alpha = .2) +
  facet_wrap(~chr) +
  ylab("slope deviation mexicana ancestry freq ~ elevation") +
  ggtitle("Where in the maize genome is mexicana ancestry best predicted by high elevation?")
ggsave("plots/mexicana_ancestry_in_maize_pops_slope_norotation_regression_elevation_wholegenome.png",
       height = 10, width = 14,
       units = "in", device = "png")

d %>%
  ggplot(aes(x = pos, y = simpleBeta_elev_maize, color = inv4m)) +
  geom_point(size = .1, alpha = .2) +
  facet_wrap(~chr) +
  ylab("slope mexicana ancestry freq ~ elevation") +
  ggtitle("Where in the maize genome is mexicana ancestry best predicted by high elevation?")
ggsave("plots/mexicana_ancestry_in_maize_pops_slope_simple_regression_elevation_wholegenome.png",
       height = 10, width = 14,
       units = "in", device = "png")


# how correlated are they?
d %>%
  dplyr::arrange(inv4m) %>%
  ggplot(aes(x = simpleBeta_elev_maize, y = zBeta_elev_maize, color = inv4m)) +
  geom_point(size = 0.1, alpha = 0.2)
ggsave("plots/slope_mexicana_in_maize_simple_regression_vs_zAnc_regression_elevation_wholegenome.png",
       height = 6, width = 6,
       units = "in", device = "png")
d %>% # perfectly correlated, just mean is shifted.
  dplyr::arrange(inv4m) %>%
  ggplot(aes(x = ancBeta_elev_maize, y = simpleBeta_elev_maize, color = inv4m)) +
  geom_point(size = 0.1, alpha = 0.2)
cor(d$ancBeta_elev_maize, d$simpleBeta_elev_maize)
d %>%
  dplyr::arrange(inv4m) %>%
  ggplot(aes(x = ancBeta_elev_maize, y = simpleBeta_elev_maize, color = inv4m)) +
  geom_point(size = 0.1, alpha = 0.2)
summary(d$ancBeta_elev_maize)
summary(d$zBeta_elev_maize) # should be centered at zero?
summary(d$simpleBeta_elev_maize) # will be positive for most loci b/c + association with mean ancestry and elevation 
# but the mean positive slope is slight


# look at ztz across the genome:
d %>%
  ggplot(aes(x = pos, y = ztz, color = inv4m)) +
  geom_point(size = .1, alpha = .2) +
  facet_wrap(~chr) +
  ylab("deviation from drift model (ztz)") +
  ggtitle("Where does mexicana ancestry in maize look least neutrally distributed across pops?")
ggsave("plots/mexicana_ancestry_in_maize_ztz_wholegenome.png",
       height = 10, width = 14,
       units = "in", device = "png")

# how many populations are participating in a 'peak'?
# let's start with the distribution of shared high mexicana across maize populations:
# by standardizing ancestry, and letting a 'peak' be anything > 2 s.d. above population mean:
maize_anc_s <- t(apply(maize_anc, 1, function(x) (x - zAnc_maize$alpha)/sqrt(diag(zAnc_maize$K)))) # normalize by subtracting pop mean and dividing by sqrt variance (= diag of K matrix)
mexicana_anc_s <- t(apply(mexicana_anc, 1, function(x) (x - zAnc_mexicana$alpha)/sqrt(diag(zAnc_mexicana$K))))
summary(maize_anc_s)
summary(maize_anc) # some populations hit 2 s.d. before getting to 100% maize ancestry; should only use this to look at high mexicana
summary(mexicana_anc_s)
summary(mexicana_anc) # all populations hit 2 s.d. before getting to 100% mexicana ancestry; can only use this to look at low mexicana/high maize

# number of populations sharing high minor-ancestry peak:
maize_shared <- apply(maize_anc_s, 1, function(x) sum(x >= 2))
hist(maize_shared[maize_shared > 0])
mexicana_shared <- apply(mexicana_anc_s, 1, function(x) sum(x <= -2))
hist(mexicana_shared[mexicana_shared > 0])
pops_shared <- data.frame(maize = apply(maize_anc_s >= 2 | maize_anc > .99, 1, sum),
                          mexicana = apply(mexicana_anc_s <= -2 | mexicana_anc < .01, 1, sum),
           sd = 2) %>%
  bind_rows(., data.frame(maize = apply(maize_anc_s >= 3 | maize_anc > .99, 1, sum),
            mexicana = apply(mexicana_anc_s <= -3 | mexicana_anc < .01, 1, sum),
            sd = 3)) %>%
  bind_rows(., data.frame(maize = apply(maize_anc_s >= 4 | maize_anc > .99, 1, sum),
                          mexicana = apply(mexicana_anc_s <= -4 | mexicana_anc < .01, 1, sum),
                          sd = 4)) %>%
  bind_rows(., data.frame(maize = apply(maize_anc_s >= 1 | maize_anc > .99, 1, sum),
                          mexicana = apply(mexicana_anc_s <= -1 | mexicana_anc < .01, 1, sum),
                          sd = 1))
pops_shared %>%
  tidyr::gather(., "zea", "n_pops_sharing", c("maize", "mexicana")) %>%
  #filter(., n_pops_sharing > 0) %>%
  ggplot(aes(x = 14 - n_pops_sharing, color = zea)) +
  #ggplot(aes(x = n_pops_sharing, color = zea)) +
  stat_ecdf(position = "identity") +
  xlab("number pops (of 14) WITHOUT elevated minor ancestry") +
  ylab("cumulative frequency across all sites w/ ancestry calls") +
  facet_wrap(~sd) +
  #ggtitle("# populations sharing elevated minor ancestry over X s.d. above pop mean") +
  #scale_x_discrete(name = "# pops sharing elevated ancestry", 
  #                   limits=c(14:0))
  ggtitle("# populations NOT sharing elevated ancestry (>99% or 1,2,3,4 s.d. above pop mean)")
ggsave("plots/cum_dist_shared_high_minor_ancestry.png",
       height = 12, width = 13, units = "in", device = "png")

pops_shared %>%
  tidyr::gather(., "zea", "n_pops_sharing", c("maize", "mexicana")) %>%
  ggplot(aes(x = n_pops_sharing, fill = zea)) +
  geom_histogram(position = "identity", alpha = .5) +
  facet_wrap(~sd) +
  ggtitle("# populations sharing elevated minor ancestry (> 99% or X s.d. above pop mean)")
ggsave("plots/hist_shared_high_minor_ancestry.png",
       height = 12, width = 13, units = "in", device = "png")

# observations: in general, maize is more likely to have no populations with > 1,2,3,4 s.d. above the mean for minor ancestry.
# this could be a power issue because in general mexicana has lower % maize than maize has % mexicana admixture.
# at the lower cutoff thresholds, maize tends to have more loci with broadly shared peaks whereas in mexicana minor ancestry peaks are 
# concentrated around 1-5 populations
# In contrast, maize has some peaks at 1 s.d. above the mean shared among all 14 maize pops, and a noticeable uptick at 8 populations


# Q: do these tend to be the same 8 populations? a cluster? Are they the highest elevation pops? Is this all the inv4m?
sort(table(apply((maize_anc_s > 2)[apply(maize_anc_s > 2, 1, sum) == 8, ], 1, paste, collapse = "-")), decreasing = T)
maize_pops[c(F, F, F, T, T, T, T, F, F, T, F, T, T, T)] # n=146 the most frequent set of 8

meta.pops[meta.pops$pop %in% maize_pops[c(F, T, T, T, F, F, T, F, T, T, T, T, F, F)], c("ELEVATION", "LOCALITY", "pop")] %>%
  arrange(ELEVATION) # n=186
meta.pops[meta.pops$pop %in% maize_pops[c(F, F, T, T, T, T, F, T, F, F, F, T, T, T)], c("ELEVATION", "LOCALITY", "pop")] %>%
  arrange(ELEVATION) # n = 390; possibly mostly the inversion
# not working -- what portion of this is the inversion?
#sort(table(apply(((maize_anc_s > 2)[apply(maize_anc_s > 2, 1, sum) == 8])[!d$inv4m, ], 
#                 1, paste, collapse = "-")), decreasing = T)


# does the pattern look different if I restrict it to just zTz outliers (presumably non-neutral sites)?



# what do these patterns of shared ancestry look like for low minor ancestry? <.01 since all pops are hitting lower bound threshold before 2 s.d.
data.frame(maize = apply(maize_anc, 1, function(x) sum(x < .01)),
           mexicana = apply(mexicana_anc, 1, function(x) sum(x > .99))) %>%
  tidyr::gather(., "zea", "n_pops_sharing") %>%
  ggplot(., aes(x = n_pops_sharing, fill = zea)) +
  geom_histogram(position = "identity", alpha = .5) +
  ylab("number of loci") +
  ggtitle("# populations sharing low (<.01) minor ancestry")
ggsave("plots/hist_shared_low_minor_ancestry_01.png",
       height = 6, width = 8, units = "in", device = "png")

data.frame(maize = apply(maize_anc < 0.1 | maize_anc_s <= -2, 1, sum),
           mexicana = apply(mexicana_anc > .99 | mexicana_anc_s >= 2, 1, sum)) %>%
  tidyr::gather(., "zea", "n_pops_sharing") %>%
  ggplot(., aes(x = n_pops_sharing, fill = zea)) +
  geom_histogram(position = "identity", alpha = .5) +
  ylab("number of loci") +
  ggtitle("# populations sharing low (< 0.01 or 2 sd below mean) minor ancestry")
ggsave("plots/hist_shared_low_minor_ancestry_01_or2sd.png",
       height =6, width = 8, units = "in", device = "png")
  

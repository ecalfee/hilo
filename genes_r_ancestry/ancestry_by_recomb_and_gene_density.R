# in this script I use recombination rate
# and gene density to predict ancestry patterns across the genome
# testing the hypothesis that minor ancestry is lowest in
# genomic regions with high gene density and low recombination
library(quantreg) # for regression on quantiles (percent ancestry)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(bedr)
library(IRanges)

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
  left_join(., unique(dplyr::select(hilo, c("popN", "zea", "symp_allo", "RI_ACCESSION", "GEOCTY", "LOCALITY"))), by = c("popN")) %>%
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
lm_coding_density_alpha <- with(filter(d2, !inv4m), 
                                lm(alpha ~ coding_density*group_alpha))
lm_coding_density_alpha_maize <- with(filter(d, !inv4m), 
                                lm(meanAlpha_maize ~ coding_density))
lm_coding_density_alpha_mex <- with(filter(d, !inv4m), 
                                      lm(meanAlpha_mex ~ coding_density))

lm_alpha <- with(filter(d2, !inv4m),
                 lm(alpha ~ group_alpha))
lm_coding_bp_alpha_maize <- with(filter(d, !inv4m),
                           lm(meanAlpha_maize ~ coding_bp))
lm_coding_bp_alpha_mex <- with(filter(d, !inv4m),
                                 lm(meanAlpha_mex ~ coding_bp))

summary(lm_coding_density_alpha)
summary(lm_alpha)
summary(lm_coding_density_alpha_maize) # highly sig. but small
summary(lm_coding_density_alpha_mex) # explains >4x the variance in ancestry as coding density does in maize
summary(lm_coding_bp_alpha_maize)
summary(lm_coding_bp_alpha_mex)

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
lm_maize <- filter(d, !inv4m) %>% # test with and without inversion
  mutate(noncoding_bp = width_bp - coding_bp) %>%
  lm(meanAlpha_maize ~ coding_bp + noncoding_bp, data = .)
lm_mex <- filter(d, !inv4m) %>%
  mutate(noncoding_bp = width_bp - coding_bp) %>%
  lm(meanAlpha_mex ~ coding_bp + noncoding_bp, data = .)
summary(lm_maize)
summary(lm_mex) # unexpected negative estimate for effect of noncoding_bp
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
# not super convincing but there are so many points in low coding density side of this fig.
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
  dplyr::select(chr, pos, inv_any, meanAlpha_maize, meanAlpha_mex) %>%
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
maize_shared <- apply(maize_anc_s, 1, function(x) sum(x >= 1))
maize_shared2 <- apply(maize_anc_s, 1, function(x) sum(x >= 2))
hist(maize_shared[maize_shared > 0])
mexicana_shared <- apply(mexicana_anc_s, 1, function(x) sum(x <= -1))
mexicana_shared2 <- apply(mexicana_anc_s, 1, function(x) sum(x <= -2))
hist(mexicana_shared[mexicana_shared > 0])

pops_shared_strict <- data.frame(maize = apply(maize_anc_s >= 2 | maize_anc > .99, 1, sum),
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
pops_shared_strict %>%
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
ggsave("plots/cum_dist_shared_high_minor_ancestry_reverse.png",
       height = 12, width = 13, units = "in", device = "png")

pops_shared_strict %>%
  tidyr::gather(., "zea", "n_pops_sharing", c("maize", "mexicana")) %>%
  #filter(., n_pops_sharing > 0) %>%
  ggplot(aes(x = n_pops_sharing, color = zea)) +
  stat_ecdf(position = "identity") +
  xlab("number pops (of 14) WITH elevated minor ancestry") +
  ylab("cumulative frequency across all sites w/ ancestry calls") +
  facet_wrap(~sd) +
  #scale_x_discrete(name = "# pops sharing elevated ancestry", 
  #                   limits=c(14:0)) +
  ggtitle("# populations sharing elevated ancestry (>99% or 1,2,3,4 s.d. above pop mean)")
ggsave("plots/cum_dist_shared_high_minor_ancestry.png",
       height = 12, width = 13, units = "in", device = "png")

pops_shared_strict %>%
  tidyr::gather(., "zea", "n_pops_sharing", c("maize", "mexicana")) %>%
  filter(., n_pops_sharing > 0) %>% # filter out loci with no population outliers
  ggplot(aes(x = n_pops_sharing, fill = zea)) +
  geom_histogram(position = "identity", alpha = .5) +
  facet_wrap(~sd) +
  ggtitle("# populations sharing elevated minor ancestry (> 99% or X s.d. above pop mean)")
ggsave("plots/hist_shared_high_minor_ancestry.png",
       height = 12, width = 13, units = "in", device = "png")

# what if I use a stringent cutoff to identify to selected loci, and then ask how many populations have > 1 s.d. elevated ancestry?
pops_shared_lenient <- do.call(rbind, 
                       lapply(1:4, function(t) 
                         bind_rows(data.frame(zea = "maize",
                                              sd = t,
                                              n_pops_sharing = maize_shared[apply(maize_anc_s >= t | maize_anc >= 0.99, 1, any)],
                                              stringsAsFactors = F),
                                   data.frame(zea = "mexicana",
                                              sd = t, 
                                              n_pops_sharing = mexicana_shared[apply(mexicana_anc_s <= -t | mexicana_anc <= 0.01, 1, any)],
                                              stringsAsFactors = F))))
pops_shared_lenient2 <- do.call(rbind, 
                               lapply(1:4, function(t) 
                                 bind_rows(data.frame(zea = "maize",
                                                      sd = t,
                                                      n_pops_sharing = maize_shared2[apply(maize_anc_s >= t | maize_anc >= 0.99, 1, any)],
                                                      stringsAsFactors = F),
                                           data.frame(zea = "mexicana",
                                                      sd = t, 
                                                      n_pops_sharing = mexicana_shared2[apply(mexicana_anc_s <= -t | mexicana_anc <= 0.01, 1, any)],
                                                      stringsAsFactors = F))))
pops_shared_lenient %>%
  filter(., sd >=2) %>%
  ggplot(aes(x = n_pops_sharing, fill = zea)) +
  geom_histogram(position = "identity", alpha = .5) +
  facet_wrap(~sd) +
  ggtitle("# populations w/ >1sd minor ancestry | 1 pop high > X s.d or 99%")
ggsave("plots/hist_shared_high_minor_ancestry_1sd_sharing.png",
       height = 12, width = 13, units = "in", device = "png")
pops_shared_lenient2 %>%
  filter(., sd >=3) %>%
  ggplot(aes(x = n_pops_sharing, fill = zea)) +
  geom_histogram(position = "identity", alpha = .5) +
  facet_wrap(~sd) +
  ggtitle("# populations w/ >2sd minor ancestry | 1 pop high > X s.d or 99%")
ggsave("plots/hist_shared_high_minor_ancestry_2sd_sharing.png",
       height = 12, width = 13, units = "in", device = "png")


# observations: in general, maize is more likely to have no populations with > 1,2,3,4 s.d. above the mean for minor ancestry.
# this could be a power issue because in general mexicana has lower % maize than maize has % mexicana admixture.
# at the lower cutoff thresholds, maize tends to have more loci with broadly shared peaks whereas in mexicana minor ancestry peaks are 
# concentrated around 1-5 populations
# In contrast, maize has some peaks at 1 s.d. above the mean shared among all 14 maize pops, and a noticeable uptick at 8 populations

# what if I limit my view to loci that appear non-neutral by zTz?
# ! problem: zTz is really zTz3_maize. I need to plot separately mexicana outliers
pops_shared %>%
  tidyr::gather(., "zea", "n_pops_sharing", c("maize", "mexicana")) %>%
  filter(., (zea == "maize" & rep(zTz3 >= quantile(zTz3, .99), 8)) |
                        (zea == "mexicana" & rep(zTz3_mex >= quantile(zTz3_mex, .99), 8))) %>%
  ggplot(aes(x = n_pops_sharing, color = zea)) +
  stat_ecdf(position = "identity") +
  xlab("number pops (of 14) WITH elevated minor ancestry") +
  ylab("cumulative frequency across all sites w/ ancestry calls") +
  facet_wrap(~sd) +
  ggtitle("top 1% zTz outliers: # populations sharing elevated ancestry (>99% or 1,2,3,4 s.d. above pop mean)")
ggsave("plots/cum_dist_shared_high_minor_ancestry_zTz_outliers99.png",
       height = 12, width = 13, units = "in", device = "png")


pops_shared %>%
  tidyr::gather(., "zea", "n_pops_sharing", c("maize", "mexicana")) %>%
  filter(., (zea == "maize" & rep(zTz3 >= quantile(zTz3, .99), 8)) |
           (zea == "mexicana" & rep(zTz3_mex >= quantile(zTz3_mex, .99), 8))) %>%
  ggplot(aes(x = n_pops_sharing, fill = zea)) +
  geom_histogram(position = "identity", alpha = .5) +
  facet_wrap(~sd) +
  ggtitle("top 1% zTz outliers: # populations sharing elevated minor ancestry (> 99% or X s.d. above pop mean)")
ggsave("plots/hist_shared_high_minor_ancestry_ztz_outliers99.png",
       height = 12, width = 13, units = "in", device = "png")

# top 10%
pops_shared %>%
  tidyr::gather(., "zea", "n_pops_sharing", c("maize", "mexicana")) %>%
  filter(., (zea == "maize" & rep(zTz3 >= quantile(zTz3, .90), 8)) |
           (zea == "mexicana" & rep(zTz3_mex >= quantile(zTz3_mex, .90), 8))) %>%
  ggplot(aes(x = n_pops_sharing, color = zea)) +
  stat_ecdf(position = "identity") +
  xlab("number pops (of 14) WITH elevated minor ancestry") +
  ylab("cumulative frequency across all sites w/ ancestry calls") +
  facet_wrap(~sd) +
  ggtitle("top 10% zTz outliers: # populations sharing elevated ancestry (>99% or 1,2,3,4 s.d. above pop mean)")
ggsave("plots/cum_dist_shared_high_minor_ancestry_zTz_outliers90.png",
       height = 12, width = 13, units = "in", device = "png")


pops_shared %>%
  tidyr::gather(., "zea", "n_pops_sharing", c("maize", "mexicana")) %>%
  filter(., (zea == "maize" & rep(zTz3 >= quantile(zTz3, .90), 8)) |
           (zea == "mexicana" & rep(zTz3_mex >= quantile(zTz3_mex, .90), 8))) %>%
  ggplot(aes(x = n_pops_sharing, fill = zea)) +
  geom_histogram(position = "identity", alpha = .5) +
  facet_wrap(~sd) +
  ggtitle("top 10% zTz outliers: # populations sharing elevated minor ancestry (> 99% or X s.d. above pop mean)")
ggsave("plots/hist_shared_high_minor_ancestry_ztz_outliers90.png",
       height = 12, width = 13, units = "in", device = "png")




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


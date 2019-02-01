# in this script I use recombination rate
# and gene density to predict ancestry patterns across the genome
# testing the hypothesis that minor ancestry is lowest in
# genomic regions with high gene density and low recombination
library(quantreg) # for regression on quantiles (percent ancestry)

# first, load data:
source("ZAnc_statistic.R") 
d = data.frame(pos, ZAnc_maize = maize[["ZAnc"]], ZAnc_mex = mex[["ZAnc"]], ZAnc_all = all[["ZAnc"]])

# first a basic plot of estimated mexicana ancestry
# by high vs. low recomb. rate regions
png(paste("../plots/alpha_by_recomb_rate.png"), # saves plot as pin in ../plots/
    height = 5, width = 7, units = "in", res = 150)
plot(y = low_r$alpha, x = high_r$alpha,
     main = "mean pop. mexicana ancestry by recomb. rate",
     ylab = "alpha low r",
     xlab = "alpha high r",
     col = ifelse(pops$zea == "mexicana", "blue", "orange"))
legend(x = "topleft",
       legend = c("maize", "mexicana", "1-to-1 line"), 
       col = c("orange", "blue", "black"),
       pch = c(1, 1, 5))
abline(a = 0, b = 1)
dev.off()

# plot ZAnc against recombination rates in 0.1cM windows around SNPs
d %>%
  filter(., chr == "4") %>%
  ggplot(., aes(x = pos, y = rate)) + 
  geom_point(cex= .1, alpha = .5, col = "black") +
  geom_smooth(col = "blue")

# testing another way of plotting:
# this is a way to retain more information with fewer data points, to create smoothed lines faster,
# but without adding noise, but the predictor variable needs to be sorted (here, position is already sorted)
# reference: https://stats.stackexchange.com/questions/108077/scatterplot-smoothing-in-r-with-big-dataset-different-methods
n = nrow(filter(d, chr == "4"))
with(filter(d, chr == "4"), smoothScatter(rate ~ pos))
k <- min(ceiling(n/256), n/2)  # Window size
kernel <- c(dnorm(seq(0, 3, length.out=k)))
kernel <- c(kernel, rep(0, n - 2*length(kernel) + 1), rev(kernel[-1]))
kernel <- kernel / sum(kernel)
smooth_rate <- Re(convolve(filter(d, chr == "4")$rate, kernel))
subsample_j <- floor(seq(1, n, k/3))
with(filter(d, chr == "4"), lines(smooth_rate[subsample_j] ~ pos[subsample_j], col = "blue"))
with(filter(d, chr == "4"), lines(lowess(smooth_rate[subsample_j] ~ pos[subsample_j], f = .1), col = "green"))
with(filter(d, chr == "4"), lines(lowess(rate[subsample_j] ~ pos[subsample_j], f = .1), col = "red"))


d %>%
  #filter(., chr == "4") %>%
  ggplot(., aes(x = rate, y = ZAnc_maize)) + 
  geom_point(cex= .1, alpha = .5, col = "black") +
  geom_smooth(col = "blue")

# group rates by quantiles for barplots:
d$bin_rate = cut(d$rate, # include full range of depth in cut
                 breaks = quantile(d$rate, p = seq(0, 1, by = .2)),
                 right = T,
                 include.lowest = T)
d$coding_density = d$coding_bp/d$width_cM
d$bin_coding_density = cut(d$coding_density,
                           breaks = quantile(d$coding_density, p = seq(0, 1, by = .2)),
                           right = T,
                           include.lowest = T)
# population mean alpha is the mean of each population's mean (so each pop is weighted equally indep. of its sample size)
# whereas mean alpha is the mean mexicana ancestry at a locus across all individuals sampled
# populations are roughly sampled equally, so they are highly correlated, but not the same
d$pop_meanAlpha_maize = colMeans(maize_anc)
d$pop_meanAlpha_mex = colMeans(mex_anc)
N_per_pop_mex = sapply(rownames(mex_anc), function(i) pops[pops$name_SHORT == i, "N_haps"]/2)
d$meanAlpha_mex = t((N_per_pop_mex %*% mex_anc)/sum(N_per_pop_mex))[ , 1]
N_per_pop_maize = sapply(rownames(maize_anc), function(i) pops[pops$name_SHORT == i, "N_haps"]/2)
d$meanAlpha_maize = t((N_per_pop_maize %*% maize_anc)/sum(N_per_pop_maize))[ , 1]


# exclude inversion area
d_split <- d %>%
  filter(., !inv_any) %>%
  gather(., "group_alpha", "alpha", 
         c("meanAlpha_mex", "meanAlpha_maize")) %>%
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
box_rate_alpha <- 
  ggplot(d_split, aes(y = alpha, x = bin_rate,
                      fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_boxplot() +
  facet_wrap(~group_alpha) +
  labs(x = "recombination rate (0.1cM winds)", 
       y = "Mean mex. ancestry by group",
       main = "Effect of local recomb. on mean mex ancestry (alpha)")
plot(box_rate_alpha)
ggsave("../plots/box_rate_alpha_0.1cM_windows.png", 
       height = 15, width = 20, 
       units = "in", device = "png")
# same but as a violin plot to show distribution:
violin_rate_alpha <- 
  ggplot(d_split, aes(y = alpha, x = bin_rate,
                      fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_violin() +
  facet_wrap(~group_alpha) +
  labs(x = "recombination rate (0.1cM winds)", 
       y = "Mean mex. ancestry by group",
       main = "Effect of local recomb. on mean mex ancestry (alpha)")
plot(violin_rate_alpha)
ggsave("../plots/violin_rate_alpha_0.1cM_windows.png",
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
  ggplot(d_split, aes(y = alpha, x = bin_coding_density,
                      fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_boxplot() +
  facet_wrap(~group_alpha) +
  labs(x = "coding density (bp/cM)", 
       y = "Mean mex. ancestry by group",
       main = "Effect of local coding density on mean mex ancestry (alpha)")
plot(box_coding_density_alpha)
ggsave("../plots/box_coding_density_alpha_0.1cM_windows.png", 
       height = 15, width = 20, 
       units = "in", device = "png")

# now for gene density: make bar plot for mean ancestry mexicna
violin_coding_density_alpha <- 
  ggplot(d_split, aes(y = alpha, x = bin_coding_density,
                      fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_violin() +
  facet_wrap(~group_alpha) +
  labs(x = "coding density (bp/cM)", 
       y = "Mean mex. ancestry by group",
       main = "Effect of local coding density on mean mex ancestry (alpha)")
plot(violin_coding_density_alpha)
ggsave("../plots/violin_coding_density_alpha_0.1cM_windows.png", 
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
  dplyr::select(chr, pos, inv_any, ZAnc_maize, ZAnc_mex, meanAlpha_maize, meanAlpha_mex) %>%
  mutate(start = floor(pos/l)*l) %>%
  group_by(chr, start) %>% # any SNPs within inversion, bin is within inversion
  summarize(inv_any = Reduce("|", inv_any), 
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
ggsave("../plots/violin_rate_alpha_10kb_windows.png",
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
ggsave("../plots/boxplot_rate_alpha_10kb_windows.png",
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
ggsave("../plots/boxplot_perc_coding_alpha_10kb_windows.png",
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
ggsave("../plots/boxplot_coding_density_alpha_10kb_windows.png",
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
ggsave("../plots/boxplot_rate_perc_coding_10kb_windows.png",
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
  dplyr::select(chr, pos, inv_any, ZAnc_maize, ZAnc_mex, meanAlpha_maize, meanAlpha_mex) %>%
  mutate(start_100kb = floor(pos/100000)*100000) %>%
  group_by(chr, start_100kb) %>% # any SNPs within inversion, bin is within inversion
  summarize(inv_any = Reduce("|", inv_any), 
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
ggsave("../plots/scatter_perc_coding_by_r_100kb_windows.png",
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
ggsave("../plots/boxplot_rate_alpha_100kb_windows.png",
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
ggsave("../plots/boxplot_perc_coding_alpha_100kb_windows.png",
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
ggsave("../plots/boxplot_coding_density_alpha_100kb_windows.png",
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

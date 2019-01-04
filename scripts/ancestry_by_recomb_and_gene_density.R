# in this script I use recombination rate
# and gene density to predict ancestry patterns across the genome
# testing the hypothesis that minor ancestry is lowest in
# genomic regions with high gene density and low recombination

# first, load data:
source("ZAnc_statistic.R") 
d = data.frame(pos, ZAnc_maize = maize[["ZAnc"]], ZAnc_mex = mex[["ZAnc"]], ZAnc_all = all[["ZAnc"]])

# first a basic plot eof stimated mexicana ancestry
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
d$meanAlpha_maize = colMeans(maize_anc)
d$meanAlpha_mex = colMeans(mex_anc)
# exclude inversion area
d$excl_inv <- d$chr == inv[inv$ID=="Inv4m", "chrom"] & 
  d$pos > inv[inv$ID == "Inv4m", "start"] - excl_inv_buffer &
  d$pos < inv[inv$ID == "Inv4m", "end"] + excl_inv_buffer
d_split <- d %>%
  filter(., !excl_inv) %>%
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

# now for gene density: make bar plot for mean ancestry mexicna
box_coding_density_alpha <- 
  ggplot(d_split, aes(y = alpha, x = bin_coding_density,
                      fill = group_alpha)) + 
  scale_fill_manual(values = colors_maize2mex[c(1,4)]) +
  geom_boxplot() +
  facet_wrap(~group_alpha) +
  labs(x = "coding density (bp/cM)", 
       y = "Mean mex. ancestry by group",
       main = "Effect of local coding desnity on mean mex ancestry (alpha)")
plot(box_coding_density_alpha)
ggsave("../plots/box_coding_density_alpha_0.1cM_windows.png", 
       height = 15, width = 20, 
       units = "in", device = "png")


# plot ZAnc against gene density in 0.1cM windows around SNPs
d %>%
  #filter(., chr == "4") %>%
  ggplot(., aes(x = perc_coding, y = ZAnc_maize)) + 
  geom_point(cex= .1, alpha = .5, col = "black") +
  geom_smooth(col = "blue")


# plot recomb. rate against gene density in 0.1cM windows around SNPs
d %>%
  #filter(., chr == "4") %>%
  ggplot(., aes(x = rate, y = perc_coding)) + 
  geom_point(cex= .1, alpha = .5, col = "black") +
  geom_smooth(col = "blue")

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

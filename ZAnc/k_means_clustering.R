# K - means clustering:

# this script uses machine learning to try to detect dominant 'clusters' in zAnc patterns in the data
# to load data needed, first run the first part of zanc_selection_models.R

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



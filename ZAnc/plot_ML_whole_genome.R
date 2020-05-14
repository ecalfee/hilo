# results of ML model comparison:
# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
source("../colors.R")

# load data 
PREFIX="pass2_alloMAIZE"

# metadata for each individual, including pop association
pop_elev <- read.table("../data/riplasm/gps_and_elevation_for_sample_sites.txt",
                       stringsAsFactors = F, header = T, sep = "\t") %>%
  dplyr::select(., popN, ELEVATION, LAT, LON)
hilo <- read.table("../samples/hilo_meta.txt", stringsAsFactors = F, header = T, sep = "\t")
# admixture proportions per pop
meta.pops <- read.table(paste0("../global_ancestry/results/NGSAdmix/", PREFIX, "/globalAdmixtureIncludedByPopN.txt"),
                        stringsAsFactors = F, header = T) %>%
  left_join(., pop_elev, by = c("popN")) %>%
  left_join(., unique(dplyr::select(hilo, c("popN", "zea", "symp_allo", "RI_ACCESSION", "GEOCTY", "LOCALITY"))), by = c("popN")) %>%
  mutate(pop = paste0("pop", popN))

# genomic sites information
# SNPs with ancestry calls
dir_sites = paste0("../local_ancestry/results/thinnedSNPs/", PREFIX)
sites <- do.call(rbind, 
                 lapply(1:10, function(i)
                   read.table(paste0(dir_sites, "/chr", i, ".var.sites"),
                              header = F, stringsAsFactors = F)))
colnames(sites) <- c("chr", "pos", "major", "minor")

# get inversions:
inv = read.table("../data/refMaize/inversions/knownInv_v4_coord.txt",
                 stringsAsFactors = F, header = T)
colnames(inv) = c("ID", "chr", "start", "end", "length")
excl_inv_buffer = 1000000 # exclude 1 megabase around inversion
inv$excl_start = inv$start - excl_inv_buffer
inv$excl_end = inv$end + excl_inv_buffer
# which sites are within the inversion?
sites <- mutate(sites, inv4m_1Mb_buffer = (chr == inv$chr[inv$ID=="inv4m"] & 
                                             pos >= inv$excl_start[inv$ID=="inv4m"] & 
                                             pos <= inv$excl_end[inv$ID == "inv4m"]),
                inv4m = (chr == inv$chr[inv$ID=="inv4m"] & 
                           pos >= inv$start[inv$ID=="inv4m"] & 
                           pos <= inv$end[inv$ID == "inv4m"])) 

# read in ML results:
null <- read.table("results/models/pass2_alloMAIZE/maize/ML_null.txt",
                   header = T, sep = "\t") %>%
  mutate(model = "null")
allpop <- read.table("results/models/pass2_alloMAIZE/maize/ML_allpop.txt",
                   header = T, sep = "\t") %>%
  mutate(model = "allpop")
elevation <- read.table("results/models/pass2_alloMAIZE/maize/ML_elevation.txt",
                  header = T, sep = "\t") %>%
  mutate(model = "elevation")
onepop <- #lapply(1:length(alpha), function(i)
  lapply(1:2, function(i)
  read.table(paste0("results/models/pass2_alloMAIZE/maize/ML_onepop", i, ".txt"),
             header = T, sep = "\t") %>%
    mutate(modelgroup = "onepop",
           pop = i,
           model = paste(modelgroup, pop)))

# load all model results (takes a while)
model_names <- c("null", "allpop", "elevation", paste0("onepop", 1:14))
models <- lapply(1:length(model_names), function(i)
  read.table(paste0("results/models/pass2_alloMAIZE/maize/ML_", 
                    model_names[i], ".txt"),
             header = T, sep = "\t") %>%
    mutate(model = model_names[i]))
model_names2 <- factor(c("Null", "All Populations", "Elevation", 
                  filter(meta.pops, zea == "maize")$LOCALITY),
                  levels = c("Null", "All Populations", "Elevation", 
                             filter(meta.pops, zea == "maize") %>%
                               arrange(., ELEVATION) %>% 
                               .$LOCALITY))


# I am unsure how to do model weights
# but I can start with the 'winning' model
aics <- bind_cols(lapply(models, function(x) x$AICc)) %>%
  data.table::setnames(model_names) 
delta_aics <- t(apply(aics, 1, function(a) a - min(a)))
weight_aics <- t(apply(delta_aics, 1, function(a) exp(-0.5*a)/sum(exp(-0.5*a))))
# model weight across the genome
weights_genome <- apply(weight_aics, 2, sum)/nrow(weight_aics)
# how many times does that model win?
wins_genome <- apply(delta_aics, 2, function(x) sum(x == 0))/nrow(weight_aics)
  # which is the winning model? # and what is the winning model AIC?
weights_genome
wins_genome
models_df <- do.call(bind_rows, lapply(1:length(models),
                                       function(i)
                                         bind_cols(sites, models[[i]]) %>%
                                         mutate(model = model_names[i],
                                                model_name = model_names2[i])))
win_models_df <- models_df %>%
  arrange(chr, pos, AICc) %>%
  filter(!duplicated(paste0(chr, pos))) # keep only lowest AIC for each locus

p_model_wins <- win_models_df %>%
  mutate(effect = ifelse(model == "null", NA, 
                         ifelse(model == "elevation",
                                ifelse(b1 < 0, "-", "+"),
                                ifelse(b < 0, "-", "+")))) %>%
  ggplot(., aes(x = model_name, fill = effect)) +
  geom_bar(aes(y = ..count../sum(..count..))) +
  theme_classic() +
  ylab("Frequency") +
  xlab("Model with lowest AICc") +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(fill = "Direction\nof effect")
ggsave("../../hilo_manuscript/figures/model_wins_barplot.png",
       plot = p_model_wins,
       height = 4, width = 5.2, 
       units = "in", dpi = 300, 
       device = "png")

ggplot(., aes(x = model, y = frequency)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
#pivot_longer(data = ., values_to = "frequency", names_to = "model")


data.frame(model = names(weights_genome), frequency = weights_genome) %>%
  ggplot(., aes(x = model, y = frequency)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
  #pivot_longer(data = ., values_to = "frequency", names_to = "model")

models[[which(model_names == "elevation")]] %>%
  ggplot(., aes(x = b0, y = b1)) + 
  geom_point()
with(models[[which(model_names == "elevation")]],
     cor(b0, b1))


# exclude single pop models:
delta_aics3 <- t(apply(aics[ , c("null", "allpop", "elevation")], 1, function(a) a - min(a)))

weight_aics3 <- t(apply(delta_aics3, 1, function(a) exp(-0.5*a)/sum(exp(-0.5*a))))

apply(weight_aics3, 2, sum)/nrow(weight_aics3)
apply(delta_aics3, 2, function(x) sum(x == 0))/nrow(weight_aics3)

# what does RSS look like under the null compared to chi-squared?
hist(models[[1]]$RSS)
# observed AncLK vs. theoretical chisq dist.
p_AncLK_chisq <- ggplot(models[[1]], aes(x = RSS)) +
  geom_histogram(aes(y = stat(density)), fill = "#999999",
                 bins = 50) +
  stat_function(fun = dchisq,
                args = list(df = length(alpha)),
                col = "black") +
  theme_classic() +
  xlab("Anc_LK")
ggsave("../../hilo_manuscript/figures/AncLK_chisq.png",
       plot = p_AncLK_chisq,
       height = 4, width = 5.2, 
       units = "in", dpi = 300, 
       device = "png")
       #device = "tiff", compression = "lzw", type = "cairo")
# null model:
p_AncLK_chisq_null <- rbind(sel_fits[[1]], 
      sel_trunc_fits[[1]]) %>%
  filter(generating_model == "null") %>%
  ggplot(., aes(x = RSS)) +
  geom_histogram(aes(y = stat(density), 
                     fill = truncated),
                 alpha = 0.6, 
                 position = "identity",
                 bins = 25) +
  stat_function(fun = dchisq,
                args = list(df = length(alpha)),
                col = "black") +
  theme_classic() +
  xlab("Anc_LK") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9"), 
                    name = element_blank(),
                    breaks = c(T, F),
                    labels = c("Truncated MVN [0, 1]", 
                               "MVN (-Inf, Inf)"))
ggsave("../../hilo_manuscript/figures/AncLK_chisq_null.png",
       plot = p_AncLK_chisq_null,
       height = 4, width = 5.2, 
       units = "in", dpi = 300, 
       device = "png")
       #device = "tiff", compression = "lzw", type = "cairo")

# MVN simulation AncLK vs. theoretical chisq dist.
p_AncLK_chisq_sim <- rbind(sel_fits[[1]], 
                           sel_trunc_fits[[1]]) %>%
  ggplot(., aes(x = RSS)) +
  geom_histogram(aes(y = stat(density), 
                     fill = truncated),
                     alpha = 0.6, 
                     position = "identity",
                     bins = 100) +
  stat_function(fun = dchisq,
                args = list(df = length(alpha)),
                col = "black") +
  facet_grid(slope ~ intercept) +
  #xlim(c(0, 50)) +
  theme_classic() +
  xlab("Anc_LK") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9"), 
                    name = element_blank(),
                    breaks = c(T, F),
                    labels = c("Truncated MVN [0, 1]", 
                             "MVN (-Inf, Inf)"))
ggsave("../../hilo_manuscript/figures/AncLK_chisq_sim.png",
       plot = p_AncLK_chisq_sim,
       height = 7, width = 7.5, 
       units = "in", dpi = 300,
       device = "png")
       #device = "tiff", compression = "lzw", type = "cairo")

# truncation simulation AncLK vs. theoretical chisq dist.
p_AncLK_chisq_sim_zoom <- rbind(sel_fits[[1]], 
                           sel_trunc_fits[[1]]) %>%
  ggplot(., aes(x = RSS)) +
  geom_histogram(aes(y = stat(density), 
                     fill = truncated),
                 alpha = 0.6, 
                 position = "identity",
                 bins = 50) +
  stat_function(fun = dchisq,
                args = list(df = length(alpha)),
                col = "black") +
  facet_grid(slope ~ intercept) +
  xlim(c(0, 50)) +
  theme_classic() +
  xlab("Anc_LK") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9"), 
                    name = element_blank(),
                    breaks = c(T, F),
                    labels = c("Truncated MVN [0, 1]", 
                               "MVN (-Inf, Inf)"))
ggsave("../../hilo_manuscript/figures/AncLK_chisq_sim_zoom.png",
       plot = p_AncLK_chisq_sim_zoom,
       height = 7, width = 7.5, 
       units = "in", dpi = 300,
       device = "png")
       #device = "tiff", compression = "lzw", type = "cairo")







# compare to ML MVN elevation model to simple linear regression:
simple_bElev <- read.table("results/models/pass2_alloMAIZE/maize/simple_bElev_anc.txt",
                           header = T, sep = "\t")
# maybe the more relevant comparison is to simple logistic not linear regression
r <- cbind(sites, elevation, simple_bElev)
r_small <- r[c(T, rep(F, 100)),]
ggplot(r_small, aes(color = inv4m, x = envWeights*1000, y = b1)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("bElev - simple linear model") +
  ylab("bElev - MVN")
ggsave("plots/bElev_MVN_vs_lm.png", height = 5, width = 5, units = "in", dpi = 300)
ggplot(r_small, aes(color = inv4m, x = X.Intercept., 
                    y = b0)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0)
cor(r_small$bElev, r_small$envWeights) # correlation = 0.86
cor(r_small$b0, r_small$X.Intercept.) # corr = 0.05
ggplot(r_small, aes(color = inv4m, x = pos, y = b1)) +
  geom_point(size = 0.1) +
  facet_wrap(~chr)
ggplot(r_small, aes(color = inv4m, x = pos, y = envWeights*1000)) +
  geom_point(size = 0.1) +
  facet_wrap(~chr)

# test out likelihood ratio test on simulated null data
# we expect:
par(mfrow=c(2,2))
#-2ln(null) + 2ln(allpop) ~ chisq_1
hist(-2*sel_fits[["null"]][sel_fits[["null"]]$generating_model == "null", "ll"] + 
       2*sel_fits[["allpop"]][sel_fits[["null"]]$generating_model == "null", "ll"], 
     freq = F, breaks = 100, main = "likelihood ratio allpop/null", xlab = "-2(diff ll)")
curve(dchisq(x, df = 1), from = 0, to = 10, lwd = 2, col = "blue", add = T)

#-2ln(null) + 2ln(1pop) ~ chisq_1
hist(-2*sel_fits[["null"]][sel_fits[["null"]]$generating_model == "null", "ll"] + 
       2*sel_fits[["onepop7"]][sel_fits[["null"]]$generating_model == "null", "ll"], 
     freq = F, breaks = 100, main = "likelihood ratio onepop7/null", xlab = "-2(diff ll)")
curve(dchisq(x, df = 1), from = 0, to = 10, lwd = 2, col = "blue", add = T)

#-2ln(null) + 2ln(elev) ~ chisq_2
hist(-2*sel_fits[["null"]][sel_fits[["null"]]$generating_model == "null", "ll"] + 
       2*sel_fits[["elevation"]][sel_fits[["null"]]$generating_model == "null", "ll"], 
     freq = F, breaks = 100, main = "likelihood ratio elev/null", xlab = "-2(diff ll)")
curve(dchisq(x, df = 2), from = 0, to = 10, lwd = 2, col = "green", add = T)

#-2ln(allpop) + 2ln(elev) ~ chisq_1
hist(-2*sel_fits[["allpop"]][sel_fits[["null"]]$generating_model == "null", "ll"] + 
       2*sel_fits[["elevation"]][sel_fits[["null"]]$generating_model == "null", "ll"], 
     freq = F, breaks = 100, main = "likelihood ratio elev/allpop", xlab = "-2(diff ll)")
curve(dchisq(x, df = 1), from = 0, to = 10, lwd = 2, col = "blue", add = T)
par(mfrow=c(1,1))


# test multiple cutoffs based on alphas for chi-sq dist against null
sig_alphas = c(0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)
bonferroni_alphas = sig_alphas/(length(models) - 1) # Bonferroni correction, don't include null in multiple-comparisons count
chisq_cv = lapply(1:2, function(k) qchisq(1 - bonferroni_alphas, df = k)) # critical values 
# calculated for 1 or 2 degrees of freedom more than null model

# make barplots for diff cutoffs decide cutoff
win_models_df2 <- bind_cols(sites, models[[1]]) %>%
  dplyr::select(chr, pos, ll) %>%
  dplyr::rename(ll_null = ll) %>%
  left_join(win_models_df, ., by = c("chr", "pos")) %>%
  mutate(lik_ratio = -2*ll_null + 2*ll) %>%
  mutate(p = pchisq(lik_ratio, lower.tail = F, df = k)) %>%
  mutate(p_bonferroni = p*(length(models) - 1)) %>% # correction for multiple testing (count all but null model)
  mutate(effect = ifelse(model == "null", NA, 
                         ifelse(model == "elevation",
                                ifelse(b1 < 0, "-", "+"),
                                ifelse(b < 0, "-", "+"))))

# filter by critical values and make a bunch of plots
# summarise data that would go into barplots
win_models_freq <- do.call(rbind, lapply(sig_alphas, function(a)
  win_models_df2 %>%
  filter(., p_bonferroni < a) %>%
  dplyr::group_by(., model, model_name, effect) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(freq = count/sum(count),
         alpha = a)))

# plot all together
p_wins_barplot_facet_alpha <- ggplot(win_models_freq, 
                                     aes(x = model_name, fill = effect)) +
  geom_bar(aes(y = freq), stat = "identity") +
  theme_classic() +
  ylab("Frequency") +
  xlab("Best model") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(fill = "Direction\nof effect") +
  facet_wrap(~alpha, nrow = 3) +
  scale_fill_manual(values = col_pos_neg)

# report what fraction of the genome is an outlier too
p_sig_alphas <- win_models_freq %>% # make log10 scale
  group_by(alpha) %>%
  summarise(count = sum(count)) %>%
  mutate(prop_sig = count/nrow(sites)) %>%
  ggplot(., aes(x = -log10(alpha), y = prop_sig * 100)) +
  geom_point() +
  ylab("% of loci sig. outliers") +
  ylim(c(0,15)) +
  xlab(expression("-log"[10]~alpha)) +
  #scale_x_reverse()+ 
  theme_light()
# grob graphs together
p_wins_alpha = grid.arrange(ncol = 2, 
             grobs = list(p_wins_barplot_facet_alpha + ggtitle("A"), 
                          p_sig_alphas + ggtitle("B")),
             widths = c(5,3))

ggsave(paste0("../../hilo_manuscript/figures/model_wins_barplot_alphas.png"),
       p_wins_alpha,
       height = 5, width = 7.5, 
       units = "in", dpi = 300, 
       device = "png")


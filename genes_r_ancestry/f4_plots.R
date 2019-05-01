library(dplyr)
library(ggplot2)
# what does f4 look like across the genome?
# distribution with recombination rate and gene density?
f4_nonadmix <- cbind(read.table("results/f4/pass2_alloMAIZE/parv_allo.maize_allo.mexicana_trip/f4_10kb.bed", stringsAsFactors = F), 
                     read.table("results/f4/pass2_alloMAIZE/parv_allo.maize_allo.mexicana_trip/SNP_count_10kb.bed", stringsAsFactors = F)$V5)
colnames(f4_nonadmix) <- c("chr", "start", "end", "10kb_region", "f4_nonadmix", "n_SNPs_nonadmix")
f4_maize <- cbind(read.table("results/f4/pass2_alloMAIZE/parv_symp.maize_allo.mexicana_trip/f4_10kb.bed", stringsAsFactors = F), 
                  read.table("results/f4/pass2_alloMAIZE/parv_symp.maize_allo.mexicana_trip/SNP_count_10kb.bed", stringsAsFactors = F)$V5)
colnames(f4_maize) <- c("chr", "start", "end", "10kb_region", "f4_maize", "n_SNPs_maize")
f4_mexicana <- cbind(read.table("results/f4/pass2_alloMAIZE/parv_symp.mexicana_allo.mexicana_trip/f4_10kb.bed", stringsAsFactors = F), 
                     read.table("results/f4/pass2_alloMAIZE/parv_symp.mexicana_allo.mexicana_trip/SNP_count_10kb.bed", stringsAsFactors = F)$V5)
colnames(f4_mexicana) <- c("chr", "start", "end", "10kb_region", "f4_mexicana", "n_SNPs_mexicana")


genome <- read.table("../data/refMaize/windows_10kb/gene_overlap.bed", stringsAsFactors = F) %>%
  cbind(., read.table("../data/refMaize/windows_10kb/whole_genome.recomb", stringsAsFactors = F))
colnames(genome) <- c("chr", "start", "end", "10kb_region", "n_cds", "length_cds", "length_10kb", "perc_coding", "r")

# combine data:
d <- left_join(genome, f4_nonadmix, by = c("chr", "start", "end", "10kb_region")) %>%
  left_join(., f4_maize, by = c("chr", "start", "end", "10kb_region")) %>%
  left_join(., f4_mexicana, by = c("chr", "start", "end", "10kb_region")) %>%
  mutate(f4_nonadmix = as.numeric(f4_nonadmix)) %>%
  mutate(f4_maize = as.numeric(f4_maize)) %>%
  mutate(f4_mexicana = as.numeric(f4_mexicana)) %>%
  mutate(alpha_mexicana = f4_mexicana/f4_nonadmix) %>% # f4-ratio to estimate % maize
  mutate(alpha_maize = f4_maize/f4_nonadmix) %>%
  mutate(coding_density = length_cds/r) %>%
  mutate(cds5 = cut(.$coding_density, breaks = c(0, .1, 1000, 10000, 100000, 1000000), right = T, include.lowest = T)) %>%
  #mutate(cds5 = cut(.$perc_coding, breaks = quantile(x = .$perc_coding, probs = seq(0, 1, by = .5)), right = T, include.lowest = T)) %>%
  mutate(r10 = cut(.$r, breaks = quantile(x = .$r, probs = seq(0, 1, by = .1)), right = T, include.lowest = T)) %>% # bin recombination rate
  mutate(r5 = cut(.$r, breaks = quantile(x = .$r, probs = seq(0, 1, by = .2)), right = T, include.lowest = T)) # 5 quantiles
# plot
hist(d$f4_maize)
hist(d$f4_mexicana)
hist(d$f4_nonadmix)
hist(d$alpha_maize)
hist(d$alpha_mexicana)
summary(d$alpha_maize)

lma <- with(d, lm(f4_maize ~ f4_nonadmix))
summary(lma)
# an overly simplistic linear model. but slope isn't 1:1 for f4_nonadmix and f4_admix 
# these are bins which is probably bad. but higher r suggests lower minor ancestry (not what we'd predict)
lmr_maize <- with(d, lm(f4_maize ~ f4_nonadmix + r + perc_coding))
summary(lmr_maize)
lmr_mex <- with(d, lm(f4_mexicana ~ f4_nonadmix + r + perc_coding))
summary(lmr_mex)
# weak model
d1 <- mutate(d, coding_density = as.numeric(length_cds)/r)
summary(with(d1, lm(f4_maize ~ f4_nonadmix + coding_density)))

d %>%
  group_by(r5) %>%
  filter(!(is.na(f4_maize) | is.na(f4_nonadmix))) %>%
  summarise(mean(f4_maize)/mean(f4_nonadmix))

d %>%
  group_by(r5) %>%
  filter(!(is.na(f4_mexicana) | is.na(f4_nonadmix))) %>%
  summarise(mean(f4_mexicana)/mean(f4_nonadmix))

mean(d[complete.cases(d[, c("f4_maize", "f4_nonadmix")]), "f4_maize"])/mean(d[complete.cases(d[, c("f4_maize", "f4_nonadmix")]), "f4_nonadmix"])
mean(d[complete.cases(d[, c("f4_mexicana", "f4_nonadmix")]), "f4_mexicana"])/mean(d[complete.cases(d[, c("f4_mexicana", "f4_nonadmix")]), "f4_nonadmix"])

mean(d$f4_maize, na.rm = T)

sum(is.na(d$f4_nonadmix))/nrow(d) # 30% of 10kb windows have no data
sum(is.na(d$f4_maize))/nrow(d) # 30% of 10kb windows have no data
sum(is.na(d$f4_mexicana))/nrow(d) # 30% of 10kb windows have no data

sum(is.na(d$f4_nonadmix) & is.na(d$f4_maize))/nrow(d) # 30% of 10kb windows have no data. So nearly 100% overlap.
sum(is.na(d$f4_nonadmix) & is.na(d$f4_mexicana))/nrow(d) # 30% of 10kb windows have no data. So nearly 100% overlap.

mean(d$f4_maize, na.rm = T)/mean(d$f4_nonadmix, na.rm = T)

# this is very different from the all.summary file results. Why would taking a mean over I'll check by region. 
maize_regions <- read.table("results/f4/pass2_alloMAIZE/parv_symp.maize_allo.mexicana_trip/all.regions", stringsAsFactors = F, header = T)

nonadmix_regions <- read.table("results/f4/pass2_alloMAIZE/parv_allo.maize_allo.mexicana_trip/all.regions", stringsAsFactors = F, header = T)
mexicana_regions<- read.table("results/f4/pass2_alloMAIZE/parv_symp.mexicana_allo.mexicana_trip/all.regions", stringsAsFactors = F, header = T)

(sum(maize_regions$f4_sum)/sum(maize_regions$n))/(sum(nonadmix_regions$f4_sum)/sum(nonadmix_regions$n)) # reasonable estimate of maize-like ancestry in maize
(sum(mexicana_regions$f4_sum)/sum(mexicana_regions$n))/(sum(nonadmix_regions$f4_sum)/sum(nonadmix_regions$n)) # also pretty reasonable

# individual regions aren't well estimated, but mean of regions is close to the genomewide mean
summary((maize_regions$f4_sum/maize_regions$n)/(nonadmix_regions$f4_sum/nonadmix_regions$n)) 
summary((mexicana_regions$f4_sum/mexicana_regions$n)/(nonadmix_regions$f4_sum/nonadmix_regions$n))

# use data for # SNPs:
d %>%
  filter(complete.cases(dplyr::select(., c("f4_maize", "f4_nonadmix")))) %>%
  summarise((sum(f4_maize*n_SNPs_maize)/sum(n_SNPs_maize))/(sum(f4_nonadmix*n_SNPs_nonadmix)/sum(n_SNPs_nonadmix))) # good, gives me genomewide result (sanity check)
# now grouping by recombination rate:
# the lowest recombination rate 10kb regions appear to have less maize-like ancestry in maize. (!) opposite what we'd expect.
d %>%
  filter(complete.cases(dplyr::select(., c("f4_maize", "f4_nonadmix")))) %>%
  group_by(r5) %>%
  summarise((sum(f4_maize*n_SNPs_maize)/sum(n_SNPs_maize))/(sum(f4_nonadmix*n_SNPs_nonadmix)/sum(n_SNPs_nonadmix))) 
# more bins, same general trend
d %>%
  filter(complete.cases(dplyr::select(., c("f4_maize", "f4_nonadmix")))) %>%
  group_by(r10) %>%
  summarise((sum(f4_maize*n_SNPs_maize)/sum(n_SNPs_maize))/(sum(f4_nonadmix*n_SNPs_nonadmix)/sum(n_SNPs_nonadmix)))
# I can't do true quantiles for coding density, but I can visualize:
d %>% # not a super clear trend
  filter(complete.cases(dplyr::select(., c("f4_maize", "f4_nonadmix")))) %>%
  group_by(cds5) %>%
  summarise((sum(f4_maize*n_SNPs_maize)/sum(n_SNPs_maize))/(sum(f4_nonadmix*n_SNPs_nonadmix)/sum(n_SNPs_nonadmix))) 


# and mexicana; pattern is not clear:
d %>%
  filter(complete.cases(dplyr::select(., c("f4_mexicana", "f4_nonadmix")))) %>%
  group_by(r5) %>%
  summarise((sum(f4_mexicana*n_SNPs_mexicana)/sum(n_SNPs_mexicana))/(sum(f4_nonadmix*n_SNPs_nonadmix)/sum(n_SNPs_nonadmix)))
d %>%
  filter(complete.cases(dplyr::select(., c("f4_mexicana", "f4_nonadmix")))) %>%
  group_by(r10) %>%
  summarise((sum(f4_mexicana*n_SNPs_mexicana)/sum(n_SNPs_mexicana))/(sum(f4_nonadmix*n_SNPs_nonadmix)/sum(n_SNPs_nonadmix)))
# are my recombination bins reasonable? appears so..
table(d$r10)
table(d$r5)

# what about effects of coding bp? just 0 or more coding bp to start
# looks like more maize-like ancestry in both maize and mexicana in windows with at least some coding bp.
# this is a very crude measure
d %>%
  filter(complete.cases(dplyr::select(., c("f4_maize", "f4_nonadmix")))) %>%
  group_by(perc_coding > 0) %>%
  summarise((sum(f4_maize*n_SNPs_maize)/sum(n_SNPs_maize))/(sum(f4_nonadmix*n_SNPs_nonadmix)/sum(n_SNPs_nonadmix)))
d %>%
  filter(complete.cases(dplyr::select(., c("f4_mexicana", "f4_nonadmix")))) %>%
  group_by(perc_coding > 0) %>%
  summarise((sum(f4_mexicana*n_SNPs_mexicana)/sum(n_SNPs_mexicana))/(sum(f4_nonadmix*n_SNPs_nonadmix)/sum(n_SNPs_nonadmix)))



# actually I should group by r first and then summarise by taking the ratio
d %>%
  dplyr::select(., c("f4_maize", "f4_nonadmix")) %>%
  .[complete.cases(.),] %>%
  summarise(mean(f4_maize)/mean(f4_nonadmix))
summary(as.numeric(f4_nonadmix$f4_nonadmix))
d %>%
  dplyr::select(., c("f4_mexicana", "f4_nonadmix")) %>%
  .[complete.cases(.),] %>%
  summarise(mean(f4_mexicana)/mean(f4_nonadmix))
# genomewide this doesn't look great so far..hmm...

ggplot(d, aes(x = "r", y = "f4_maize")) +
  geom_point()


# how well do global estimates of admixture (NGSAdmix) match f4-ratio statistic results genomewide?
# read in population-level NGSadmix results:
ngsadmix <- read.table("../global_ancestry/results/NGSAdmix/pass2_alloMAIZE/globalAdmixtureIncludedByPopN.txt", stringsAsFactors = F, header = T)
ngsadmixInd <- read.table("../global_ancestry/results/NGSAdmix/pass2_alloMAIZE/globalAdmixtureByIncludedIndividual.txt", stringsAsFactors = F, header = T)
ngsadmixInd$zea <- ifelse(ngsadmixInd$popN > 100, "maize", "mexicana")
ngsadmix$zea <- ifelse(ngsadmix$popN > 100, "maize", "mexicana")

# group means aren't exactly the same but similar to genomewide f4-ratio results. 
# exact match isn't expected: f4 uses more SNPs and overall diff. approach; but reassuring the f4 isn't way lower
ngsadmixInd %>%
  group_by(zea) %>%
  summarise(mean(alpha_maize)) 

# I can also look at each individual population:
pops <- unique(ngsadmix$popN)
f4_pops <- lapply(pops, function(p) cbind(read.table(paste0("results/f4/pass2_alloMAIZE/parv_pop", p, "_allo.mexicana_trip/f4_10kb.bed"), stringsAsFactors = F), 
                                             read.table(paste0("results/f4/pass2_alloMAIZE/parv_pop", p, "_allo.mexicana_trip/SNP_count_10kb.bed"), stringsAsFactors = F)$V5) %>%
                    setNames(., c("chr", "start", "end", "10kb_region", "f4", "n_SNPs")) %>%
                    left_join(genome, ., by = c("chr", "start", "end", "10kb_region")) %>%
                    left_join(., f4_nonadmix, by = c("chr", "start", "end", "10kb_region")) %>%
                    mutate(f4 = as.numeric(f4)) %>%
                    mutate(f4_nonadmix = as.numeric(f4_nonadmix)) %>%
                    mutate(r10 = cut(.$r, breaks = quantile(x = .$r, probs = seq(0, 1, by = .1)), right = T, include.lowest = T)) %>% # bin recombination rate
                    mutate(r5 = cut(.$r, breaks = quantile(x = .$r, probs = seq(0, 1, by = .2)), right = T, include.lowest = T)) # 5 quantiles
)
f4_pop_alpha_est <- sapply(f4_pops, function(x) {
  x %>%
    filter(complete.cases(dplyr::select(., c("f4", "f4_nonadmix")))) %>%
    summarise((sum(f4*n_SNPs)/sum(n_SNPs))/(sum(f4_nonadmix*n_SNPs_nonadmix)/sum(n_SNPs_nonadmix)))
})
ngsadmix$f4_alpha_maize <- unlist(f4_pop_alpha_est)
ngsadmix %>%
  ggplot(aes(x = alpha_maize, y = f4_alpha_maize, size = n, color = zea)) +
  geom_point() +
  xlab("NGSadmix estimate of maize admixture proportion") +
  ylab("f4-ratio estimate of maize admixture proportion") +
  ggtitle("Comparing f4 and NGSadmix population admixture estimates") +
  geom_abline(intercept = 0, slope = 1, color = "black")
ggsave("plots/f4_vs_ngsadmix_estimate_of_alpha_maize_pass2_alloMAIZE.png",
       device = "png",
       width = 8, height = 6, units = "in",
       dpi = 300)

f4_pop_alpha_by_r_est <- lapply(f4_pops, function(x) {
  x %>%
    filter(complete.cases(dplyr::select(., c("f4", "f4_nonadmix")))) %>%
    group_by(r5) %>%
    summarise((sum(f4*n_SNPs)/sum(n_SNPs))/(sum(f4_nonadmix*n_SNPs_nonadmix)/sum(n_SNPs_nonadmix)))
})
for (i in 1:length(pops)){
  print(paste("pop", pops[i]))
  print(f4_pop_alpha_by_r_est[[i]][, 2])
}

f4_pop_alpha_by_r10_est <- lapply(f4_pops, function(x) {
  x %>%
    filter(complete.cases(dplyr::select(., c("f4", "f4_nonadmix")))) %>%
    group_by(r10) %>%
    summarise((sum(f4*n_SNPs)/sum(n_SNPs))/(sum(f4_nonadmix*n_SNPs_nonadmix)/sum(n_SNPs_nonadmix)))
})
for (i in 1:length(pops)){
  print(paste("pop", pops[i]))
  print(f4_pop_alpha_by_r10_est[[i]][, 2])
}
# if any pattern is here in the f4-ratio stats it looks like more minoar ancestry in the lowest recomb. bins over 10kb windows


# exclude chromosome 4 due to the inversion:
d %>% # same pattern
  filter(complete.cases(dplyr::select(., c("f4_maize", "f4_nonadmix")))) %>%
  filter(chr != 4) %>%
  group_by(r5) %>%
  summarise((sum(f4_maize*n_SNPs_maize)/sum(n_SNPs_maize))/(sum(f4_nonadmix*n_SNPs_nonadmix)/sum(n_SNPs_nonadmix))) 
d %>% # weak pattern
  filter(complete.cases(dplyr::select(., c("f4_mexicana", "f4_nonadmix")))) %>%
  filter(chr != 4) %>%
  group_by(r5) %>%
  summarise((sum(f4_mexicana*n_SNPs_mexicana)/sum(n_SNPs_mexicana))/(sum(f4_nonadmix*n_SNPs_nonadmix)/sum(n_SNPs_nonadmix))) 

# maybe these estimates are just too unstable and I should be using angsd's built-in -doAbbaBaba2 f4 estimator..

library(dplyr)
library(ggplot2)
# what does f4 look like across the genome?
# distribution with recombination rate and gene density?
f4_nonadmix <- read.table("results/f4/all/parv_allo.maize_allo.mexicana_trip/f4_10kb.bed", stringsAsFactors = F)
colnames(f4_nonadmix) <- c("chr", "start", "end", "10kb_region", "f4_nonadmix")
f4_maize <- read.table("results/f4/all/parv_symp.maize_allo.mexicana_trip/f4_10kb.bed", stringsAsFactors = F)
colnames(f4_maize) <- c("chr", "start", "end", "10kb_region", "f4_maize")
f4_mexicana <- read.table("results/f4/all/parv_symp.mexicana_allo.mexicana_trip/f4_10kb.bed", stringsAsFactors = F)
colnames(f4_mexicana) <- c("chr", "start", "end", "10kb_region", "f4_mexicana")


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
  mutate(alpha_maize = f4_maize/f4_nonadmix)
  

# plot
hist(d$f4_maize)
hist(d$f4_mexicana)
hist(d$f4_nonadmix)
hist(d$alpha_maize)
hist(d$alpha_mexicana)
summary(d$alpha_maize)


# actually I should group by r first and then summarise by taking the ratio
d %>%
  dplyr::select(., c("f4_maize", "f4_nonadmix")) %>%
  .[complete.cases(.),] %>%
  summarise(mean(f4_maize)/mean(f4_nonadmix))
d %>%
  dplyr::select(., c("f4_mexicana", "f4_nonadmix")) %>%
  .[complete.cases(.),] %>%
  summarise(mean(f4_mexicana)/mean(f4_nonadmix))
# genomewide this doesn't look great so far..hmm...

ggplot(d, aes(x = "r", y = "f4_maize")) +
  geom_point()


# how well do global estimates of admixture (NGSAdmix) match f4-ratio statistic results genomewide?
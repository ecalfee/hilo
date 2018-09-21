# this script plots posterior results from ancestry_hmm
library(dplyr)
library(tidyr)
library(ggplot2)

dir_in = "../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/"
alphas <- read.table("../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/input/globalAdmixtureByPopN.txt",
                     stringsAsFactors = F, header = F)
colnames(alphas) <- c("popN", "alpha_maize", "alpha_mex")
meta <- read.table("../data/pass1_ids.txt", stringsAsFactors = F, 
                   header = T, sep = "\t") %>%
  filter(est_coverage >= 0.05) %>%
  group_by(popN, zea, symp_allo, RI_ACCESSION, GEOCTY, LOCALITY) %>%
  summarise(tot_pop_coverage = sum(est_coverage)) %>%
  left_join(., alphas, by = c("popN"))

pops = meta[(meta$symp_allo == "sympatric" & 
              meta$alpha_maize > 0 & 
              meta$alpha_mex > 0), "popN"] %>%
  unlist(.) # won't work as tibble in lapply below

get_boots = function(path){
  boots = do.call(rbind, lapply(pops, function(pop) read.table(paste0(path, 
                                                                      "bootstrap_pop", 
                                                                      pop, ".txt"), 
                                                               stringsAsFactors = F, 
                                                               header = F, skip = 1)$V2))
  colnames(boots)<- paste0("boot", 1:11)
  boots = as.data.frame(boots)
  boots$popN = pops
  boots = gather(boots, "bootstrap", "t_est", 1:11) %>%
    left_join(., meta, by = c("popN"))
  return(boots)
}
boots10 = get_boots(path = paste0(dir_in, "output_bootT/"))
boots10$Ne = 10000
boots5 = get_boots(path = paste0(dir_in, "output_bootT_Ne5k/"))
boots5$Ne = 5000
boots50 = get_boots(path = paste0(dir_in, "output_bootT_Ne50k/"))
boots50$Ne = 50000
boots = rbind(boots10, boots5, boots50)


# no clear relationship between age of admixture estimate (t)
# and total population coverage or alpha
p1 = ggplot(boots10, aes(col = zea)) +
  geom_point(aes(x = tot_pop_coverage, y = t_est)) +
  ggtitle("time estimate by total population coverage")
plot(p1)
p1
ggsave("../plots/t_est_by_pop_coverage.png", plot = p1, device = png(), 
       width = 7, height = 5, units = "in",
       dpi = 200)
p2 = ggplot(filter(meta, symp_allo=="sympatric"), 
       aes(col = zea)) +
  geom_point(aes(x = tot_pop_coverage, y = alpha_mex)) +
  ggtitle("alpha mean admixture estimate by total population coverage")
p2
ggsave("../plots/alpha_est_by_pop_coverage.png", plot = p2, device = png(), 
       width = 7, height = 5, units = "in",
       dpi = 200)

p3 = ggplot(boots10, aes(x = zea, y = t_est, color = zea)) +
  geom_violin() + 
  geom_point(size=1, position = position_jitter(w=0.05)) +
  theme_minimal() + 
  facet_wrap(~LOCALITY) # show plot
plot(p3)
ggsave("../plots/t_est_by_location_violin.png", plot = p3, 
       device = png(), 
       width = 7, height = 5, units = "in",
       dpi = 200)

p4 = ggplot(boots10, aes(x = zea, y = t_est, color = zea)) +
  geom_point() +
  ggtitle("10 block bootstrap time of admixture estimates", 
          sub = "5 pops have no pair b/c no inferred admixture") +
  theme_minimal() + 
  facet_wrap(~LOCALITY) # show plot
ggsave("../plots/t_est_by_location_points.png", plot = p4, 
       device = png(), 
       width = 7, height = 5, units = "in",
       dpi = 200)

# now with different Ne estimates too:
p5 = ggplot(boots, aes(x = zea, y = t_est, color = as.factor(Ne))) +
  geom_point(size=1, position = position_jitter(w=0.2)) + # jitter to see overlap
  ggtitle("Vary Ne: 10 block bootstrap time of admixture estimates", 
          sub = "5 pops have no pair b/c no inferred admixture") +
  theme_minimal() + 
  facet_wrap(~LOCALITY) # show plot
ggsave("../plots/t_est_by_location_varying_Ne.png", plot = p5, 
       device = png(), 
       width = 7, height = 5, units = "in",
       dpi = 200)


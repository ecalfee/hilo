library(dplyr)
library(ggplot2)
library(scales)
# helper file with some useful functions
source("../../covAncestry/forqs_sim/k_matrix.R") # import useful functions

dir_in = "../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/"
dir_anc = paste0(dir_in, "output_noBoot/anc/")
# admixture proportions per pop
alphas <- read.table("../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/input/globalAdmixtureByPopN.txt",
                     stringsAsFactors = F, header = F)
colnames(alphas) <- c("popN", "alpha_maize", "alpha_mex")
# metadata for each individual, including pop association
meta <- read.table("../data/pass1_ids.txt", stringsAsFactors = F, 
                   header = T, sep = "\t") %>%
  filter(est_coverage >= 0.05) %>%
  left_join(., alphas, by = c("popN"))
# individual by individual K matrix
included_inds = meta %>%
  filter(alpha_maize > 0 & alpha_mex > 0)
pops = unique(included_inds[order(included_inds$popN), c("popN", "zea", "LOCALITY", "alpha_maize", "alpha_mex")])
write.table(pops$popN, paste0(dir_in, "/included_pops.txt"),
            col.names = F, row.names = F, quote = F)

# read in input files from calc_genomewide_pop_anc_freq.R
pop_anc_list = lapply(pops$popN, function(pop) read.table(paste0(dir_anc, "/pop", pop, ".anc.freq"), 
                                                stringsAsFactors = F))

# combine population ancestry frequencies into a matrix
# where rows are populations and columns are snps
all_anc = t(do.call(cbind, pop_anc_list))
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
all = make_calcs(pop_anc = all_anc)
# just maize separately:
maize_anc = pop_anc[pops$zea == "maize", ]
maize = make_calcs(maize_anc)
# just mexicana separately:
mex_anc = pop_anc[pops$zea == "mexicana", ]
mex = make_calcs(mex_anc)

# plot K as a correlation matrix
png(paste("../plots/K_corrplot_combined.png"), # saves plot as pin in ../plots/
    height = 7, width = 7, units = "in", res = 150)
corrplot(cov2cor(pop_K), method="shade", main = "combined mex/maize")
dev.off()

sum(maize[["p_ZAnc"]] < .01)/sum(maize[["ZAnc"]] > 0)
# ~ 9% of teosinte-biased ancestry segments appear to be selected for mexicana upon rough approximation
maize_bound = qnorm(p = .01, mean = 0, sd = sqrt(nrow(maize_anc)), lower.tail = F)
sum(abs(maize[["ZAnc"]]) < maize_bound)/ncol(maize_anc) # 3.6% genome overall appears to be under selection
sum(maize[["ZAnc"]] > maize_bound)/ncol(maize_anc) # nearly all of it in the teosinte direction
mex_bound = qnorm(p = .01, mean = 0, sd = sqrt(nrow(mex_anc)), lower.tail = F)
all_bound = qnorm(p = .01, mean = 0, sd = sqrt(nrow(all_anc)), lower.tail = F)
sum(mex[["ZAnc"]] > mex_bound)/ncol(mex_anc) # hitting boundary?
sum(mex[["ZAnc"]] < -mex_bound)/ncol(mex_anc) # 1-2% maize biased


# histogram of ZAnc statistic
pop_results = list(all, maize, mex)
for (i in pop_results){
  png(paste("../plots/ZAnc_hist_", deparse(substitute(i)) , ".png"), # saves plot as pin in ../plots/
      height = 5, width = 5, units = "in", res = 150)
  hist(i[["ZAnc"]], main = paste("ZAnc", deparse(substitute(i))))
  mtext(paste0("p<.01 outliers maize:", round(sum(i[["ZAnc"]] < -all_bound)/ncol(all_anc), 3), 
               " mex:", round(sum(i[["ZAnc"]] > all_bound)/ncol(all_anc), 3))) # I suspect outliers are driven mainly by maize
  abline(v = c(all_bound, -all_bound), col = "red")
  dev.off()
}


png(paste("../plots/ZAnc_hist_combined.png"), # saves plot as pin in ../plots/
    height = 5, width = 5, units = "in", res = 150)
hist(all[["ZAnc"]], main = "Combined ZAnc")
mtext(paste0("p<.01 outliers maize:", round(sum(all[["ZAnc"]] < -all_bound)/ncol(all_anc), 3), 
                  " mex:", round(sum(all[["ZAnc"]] > all_bound)/ncol(all_anc), 3))) # I suspect outliers are driven mainly by maize
abline(v = c(all_bound, -all_bound), col = "red")
dev.off()
png(paste("../plots/ZAnc_hist_maize.png"), # saves plot as pin in ../plots/
    height = 5, width = 5, units = "in", res = 150)
hist(maize[["ZAnc"]], main = "maize ZAnc")
mtext(paste0("p<.01 outliers maize:", round(sum(maize[["ZAnc"]] < -maize_bound)/ncol(maize_anc), 3), 
             " mex:", round(sum(maize[["ZAnc"]] > maize_bound)/ncol(maize_anc), 3)))
abline(v = c(maize_bound, -maize_bound), col = "red")
dev.off()
png(paste("../plots/ZAnc_hist_mex.png"), # saves plot as pin in ../plots/
    height = 5, width = 5, units = "in", res = 150)
hist(mex[["ZAnc"]], main = "mex ZAnc")
mtext(paste0("p<.01 outliers maize:", round(sum(mex[["ZAnc"]] < -mex_bound)/ncol(mex_anc), 3), 
             " mex:", round(sum(mex[["ZAnc"]] > mex_bound)/ncol(mex_anc), 3)))
abline(v = c(mex_bound, -mex_bound), col = "red") # may be too close to the boundary to get high mex. outliers
dev.off()
plot(maize[["ZAnc"]], maize[["p_ZAnc"]], cex = .1, 
     col = ifelse(maize[["p_ZAnc"]] < .01 | maize[["p_ZAnc"]] > .99, "red", "blue"),
     main = "p_values for ZAnc statistic across sites in maize")


# plot alpha NGSAdmix vs. alpha from ancestry_hmm
pops$alpha_ancestry_hmm = all[["alpha"]]
ggplot(pops, aes(x = alpha_mex, y = alpha_ancestry_hmm)) +
  geom_point(aes(col = LOCALITY)) +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  ggtitle("genome-wide mexicana ancestry: NGSadmix vs. ancestry_hmm")
ggsave("../plots/alpha_est_NGSadmix_vs_ancestry_hmm.png", height = 5, width = 5, units = "in", device = "png")

# now get positions so I can plot ZAnc by chromosome:
pos = read.table(paste0(dir_in, "/output_noBoot/HILO1.posterior"),
                   stringsAsFactors = F, header = T)[ , 1:2]
d = data.frame(pos, ZAnc_maize = maize[["ZAnc"]], ZAnc_mex = mex[["ZAnc"]], ZAnc_all = all[["ZAnc"]])
d %>%
  #filter(., chrom == "4") %>%
  ggplot(., aes(x = position, y = ZAnc_maize)) + 
  geom_point(cex= .1, alpha = .5, col = "orange") + facet_wrap(~chrom)
ggsave("../plots/ZAnc_maize_whole_genome.png", height = 15, width = 15, units = "in", device = "png")
d %>%
  #filter(., chrom == "4") %>%
  ggplot(., aes(x = position, y = ZAnc_mex)) + 
  geom_point(cex= .1, alpha = .5, col = "blue") + facet_wrap(~chrom)
ggsave("../plots/ZAnc_mex_whole_genome.png", height = 15, width = 15, units = "in", device = "png")
d %>%
  #filter(., chrom == "4") %>%
  ggplot(., aes(x = position, y = ZAnc_all)) + 
  geom_point(cex= .1, alpha = .5, col = "black") + facet_wrap(~chrom)
ggsave("../plots/ZAnc_combined_whole_genome.png", height = 15, width = 15, units = "in", device = "png")

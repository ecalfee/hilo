# this script plots posterior results from ancestry_hmm
library(dplyr)
library(tidyr)
library(ggplot2)

pass1 <- read.table("../data/pass1_ids.txt", stringsAsFactors = F, 
                                     header = T, sep = "\t")
pop_admix_global <- read.table("../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/input/globalAdmixtureByPopN.txt",
                               stringsAsFactors = F, header = F, sep = "\t")
colnames(pop_admix_global) = c("popN", "popMaize", "popMex")
dir_in = "../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/output_noBoot/"
#popN = 23
include = left_join(pass1, pop_admix_global, by = "popN") %>%
  filter(., popMaize > 0 & popMex > 0) %>% # admixed pops
  filter(., est_coverage > 0.05) %>%
  filter(., symp_allo == "sympatric")
#ids_list = read.table(paste0("../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/input/pop",
#                  popN, ".anc_hmm.ids"), stringsAsFactors = F)$V1
#ids_list = c(1, 3, 14, 74, 78, 95, 100, 109, 130) # pop366

# get posterior for individuals at some resolution
# by default returns posterior for every 100th snp
getPostSmall = function(id, dir = dir_in, nSkip = 99){ 
  post = read.table(paste0(dir_in, "/HILO", id, ".posterior"), stringsAsFactors = F, header = T) %>%
    rename(., mex.mex = X0.2) %>%
    rename(., maize.mex = X1.1) %>%
    rename(., maize.maize = X2.0) %>%
    mutate(., hilo_id = id)
  post$max_p = apply(post, 1, function(i) max(i[3:5]))
  small1 = tidyr::gather(post[c(T, rep(F, nSkip)), ], "anc", "p", 3:5) %>%
    filter(., p == max_p) # only keep highest posterior prob. ancestry
  return(small1)
}
# get full posterior but only for a small range of values
getPostSeg = function(id, chr, start, end, buffer = 0,
                      dir = dir_in, nSkip = 0){ 
  post = read.table(paste0(dir_in, "/HILO", id, ".posterior"), stringsAsFactors = F, header = T) %>%
    filter(., chrom == chr) %>%
    filter(., position >= start - buffer) %>%
    filter(., position <= end + buffer) %>%
    rename(., mex.mex = X0.2) %>%
    rename(., maize.mex = X1.1) %>%
    rename(., maize.maize = X2.0) %>%
    mutate(., hilo_id = id)
  post$max_p = apply(post, 1, function(i) max(i[3:5]))
  seg1 = tidyr::gather(post[c(T, rep(F, nSkip)), ], "anc", "p", 3:5) %>%
    filter(., p == max_p) # only keep highest posterior prob. ancestry
  return(seg1)
}

# plot one individual skipping few vs. many
getPostSmall(id = 14, nSkip = 30) %>%
  filter(chrom == 8) %>%
  ggplot(aes(x=position, y = hilo_id)) +
  geom_point(aes(color = anc, alpha = p), size = .5)
getPostSmall(id = 14, nSkip = 99) %>% # default
  filter(chrom == 8) %>%
  ggplot(aes(x=position, y = hilo_id)) +
  geom_point(aes(color = anc, alpha = p), size = .5)

# plot a whole population
pop_23 = do.call(rbind, 
                 lapply(ids_list, function(i) 
                   getPostSmall(id=i)))
pop_23 %>%
  filter(chrom == 7) %>%
  ggplot(aes(x=position, y = hilo_id)) +
  geom_point(aes(color = anc, alpha = p), size = .5)

# for just a couple peaks and the inversions
# plot all the mexicana individuals 
# and all of the maize individuals
# check for anything weird
include_maize = include %>%
  filter(zea == "maize")
include_mex = include %>%
  filter(zea == "mexicana")
# I could make this more reasonable

maize = do.call(rbind,
                lapply(include_maize$n, function(i) 
  getPostSmall(id = i)))
#write.table(maize, paste0(dir_in, "/maize_posteriors_small.txt"),
#                          sep = "\t", col.names = T, row.names = F,
#                          quote = F)
maize = read.table(paste0(dir_in, "/maize_posteriors_small.txt"),
                   sep = "\t", header = T, stringsAsFactors = F)
mex = do.call(lapply(include_mex$n, function(i)
  getPostSmall(id = i)))
#write.table(mex, paste0(dir_in, "/mex_posteriors_small.txt"),
#            sep = "\t", col.names = T, row.names = F,
#            quote = F)
mex = read.table(paste0(dir_in, "/mex_posteriors_small.txt"),
                   sep = "\t", header = T, stringsAsFactors = F)


# plot population ancestry (with uncertainty) over some region
maize %>%
  filter(chrom == 4) %>%
  ggplot(aes(x=position, y = hilo_id)) +
  geom_point(aes(color = anc, alpha = p), size = .5)

mex %>%
  filter(chrom == 4) %>%
  ggplot(aes(x=position, y = hilo_id)) +
  geom_point(aes(color = anc, alpha = p), size = .5)

rbind(maize, mex) %>%
  filter(chrom == 4) %>%
  ggplot(aes(x=position, y = hilo_id)) +
  geom_point(aes(color = anc, alpha = p), size = .5)


plot_posteriors = function(post, chr, start, 
                           end, title, save = F){
  pA = post %>%
    filter(chrom == chr) %>%
    filter(position >= start) %>%
    filter(position <= end) %>%
    ggplot(aes(x=position, y = hilo_id)) +
    geom_point(aes(color = anc, alpha = p), size = .5) +
    ggtitle(title)
  pA
  if (save){
    ggsave(filename = paste0("../plots/post_", title, ".png"),
           plot = pA, 
           device = "png", height = 8, width = 14, units = "in")
  }
  return(pA)
}
# save some plots
plot_posteriors(post = mex, chr = 4, start = 0, end = 1000000000000,
                title = "mexicana_chr_4", save = T)
plot_posteriors(post = maize, chr = 4, start = 0, end = 1000000000000,
                title = "maize_chr_4", save = T)
plot_posteriors(post = rbind(maize, mex), chr = 4, start = 0, end = 1000000000000,
                title = "mex_and_maize_chr_4", save = T)
# plot some zAnc peaks
#4:6000000-10000000
#4:23000000-28000000
plot_posteriors(post = rbind(maize), chr = 4, 
                start = 6000000, end = 10000000,
                title = "maize_chr_4_peak1", save = T)
plot_posteriors(post = rbind(maize), chr = 4, 
                start = 23000000, end = 28000000,
                title = "maize_chr_4_peak2", save = T)
# more zoomed in (fewer points skipped):
peak1_maize = do.call(rbind, 
                      lapply(include_maize$n, 
                             function(i) 
                getPostSeg(id = i, chr = 4, start = 6000000,
                           end = 10000000, buffer = 1000000,
                           nSkip = 9))) # every 10th point plot
plot_posteriors(post = peak1_maize, chr = 4, 
                start = 0, end = 1000000000,
                title = "maize_chr_4_peak1_zoom", save = T)
# encouraging -- as we'd expect (not low confidence or all hets)
# I could additionally sort these by population
peak2_maize = do.call(rbind, 
                      lapply(include_maize$n, 
                             function(i) 
                               getPostSeg(id = i, chr = 4, start = 23000000,
                                          end = 28000000, buffer = 1000000,
                                          nSkip = 9))) # every 10th point plot
plot_posteriors(post = peak2_maize, chr = 4, 
                start = 0, end = 1000000000,
                title = "maize_chr_4_peak2_zoom", save = T)

# plot all of the inversions
inv = read.table("../data/refMaize/inversions/knownInv.txt",
                 stringsAsFactors = F, header = T)

plot_inversion = function(post, chr, start, 
                           end, ID, 
                          context = 100000000, 
                          save = T){
  pI = post %>%
    filter(chrom == chr) %>%
    filter(position >= (start - context)) %>%
    filter(position <= (end + context)) %>%
    ggplot(aes(x=position, y = hilo_id)) +
    geom_point(aes(color = anc, alpha = p), size = .5) +
    ggtitle(ID) +
    geom_vline(xintercept = start) +
    geom_vline(xintercept = end)
  pI

  if (save){
    ggsave(filename = paste0("../plots/inv_post_", ID, ".png"),
           plot = pI, 
           device = "png", height = 8, width = 14, units = "in")
  }
  return(pI)
}
# save plots of all the inversions
inv_plots = lapply(1:nrow(inv), function(i)
  plot_inversion(post = rbind(maize, mex), 
                 chr = inv$chrom[i], 
                 start = inv$start[i], 
                 end = inv$end[i], 
                 ID = paste(inv$ID[i], 
                            "maize_mex", sep = "_")))



# function to plot one individual at high definition
plot_hilo_individual <- function(id, dir_in = dir_in){
  coverage = round(pass1[pass1$n == id, "est_coverage"], 2)
  post = read.table(paste0(dir_in, "/HILO", id, ".posterior"), stringsAsFactors = F, header = T) %>%
    rename(., mex.mex = X0.2) %>%
    rename(., maize.mex = X1.1) %>%
    rename(., maize.maize = X2.0)
  post1 = tidyr::gather(post, "anc", "p", 3:5)
  small1 = tidyr::gather(post[c(T, rep(F, 99)), ], "anc", "p", 3:5)
  # percent sites inferred with high confidence
  # /3 because there are 3 rows (and hmm states) 
  # per SNP and only one 
  # can possibly be over 50% posterior
  over9 = round(sum(small1$p > .9)/(nrow(small1)/3),2)
  over8 = round(sum(small1$p > .8)/(nrow(small1)/3),2)
  # estimated global ancestry from hmm
  perc_mex_est = small1 %>% 
    group_by(., anc) %>%
    summarise(., alpha = mean(p)) %>%
    .$alpha %*% (0:2)/2
  perc_mex_est = round(perc_mex_est, 2)
  
  # plot mexicana ancestry genomewide
  p1 = small1 %>%
    filter(anc == "mex.mex") %>%
    ggplot(aes(x=position, y = p)) +
    geom_point(aes(color = anc)) +
    ggtitle(paste0("posterior prob. HILO", id, " est. cov = ", coverage), 
            paste0("est: % mex=", perc_mex_est, " % sites > .9 post= ", over9))
  p1 + facet_wrap(~chrom)
  ggsave(filename = paste0("../plots/tractsMexAnc_HILO", id, "_allCHR"),
         plot = p1 + facet_wrap(~chrom), 
         device = "png", height = 12, width = 20, units = "in")
  
  # plot posterior for mexicana 1 or 2 alleles genomewide
  p2 = small1 %>%
    filter(anc != "maize.maize") %>%
    ggplot(aes(x=position, y = p)) +
    geom_point(aes(color = anc)) +
    ggtitle(paste0("posterior prob. HILO", id, " est. cov = ", coverage), 
            paste0("est: % mex=", perc_mex_est, " % sites > .9 post= ", over9))
  p2 + facet_wrap(~chrom)
  ggsave(filename = paste0("../plots/tractsMexMaizeAnc_HILO", id, "_allCHR"),
         plot = p2 + facet_wrap(~chrom), 
         device = "png", height = 12, width = 20, units = "in")
  
  p3 = post1 %>%
    filter(., chrom == 3,
           position < 1000000) %>%
    ggplot(., aes(x=position, y = p)) +
    geom_point(aes(color = anc)) + 
    ggtitle(paste0("zoomed in posterior prob. HILO", id, " est. cov = ", coverage), 
            paste0("est: % mex=", perc_mex_est, " % sites > .9 post= ", over9))
  p3
  ggsave(filename = paste0("../plots/tractsMexMaizeAnc_HILO", id, "_zoomCHR3"),
         plot = p3, 
         device = "png", height = 12, width = 20, units = "in")
    
  p4 = small1 %>%
    filter(., chrom == 4 &
           position > 1e8*1.6 &
             position < 1e8*1.9) %>%
    ggplot(., aes(x=position, y = p)) +
    geom_point(aes(color = anc)) +
    ggtitle(paste0("zoomed in posterior prob. HILO", id, " est. cov = ", coverage), 
            paste0("est: % mex=", perc_mex_est, " % sites > .9 post= ", over9))
  p4
  ggsave(filename = paste0("../plots/tractsMexMaizeAnc_HILO", id, "_zoomCHR4"),
         plot = p4, 
         device = "png", height = 12, width = 20, units = "in")
}
# run code below to plot one individual at high definition
# for (i in ids_list) plot_hilo_individual(id = i, dir_in = dir_in)



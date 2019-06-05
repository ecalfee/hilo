# this script plots posterior results from ancestry_hmm, summarised as mean ancestry per population
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

# allopatric haplotypes to compare to
allo_haps <- c("maize.allo.4Low16", "mexicana.allo.withXochi35")

# get posterior for individuals at some resolution
# by default returns posterior for every 100th snp
getPostSmall = function(id, dir = dir_in, nSkip = 99){ 
  post = read.table(paste0(dir_in, "/HILO", id, ".posterior"), stringsAsFactors = F, header = T) %>%
    rename(., mex.mex = X0.2) %>%
    rename(., maize.mex = X1.1) %>%
    rename(., maize.maize = X2.0) %>%
    mutate(., hilo_id = id)
  post$max_p = apply(post[ , 3:5], 1, max)
  small1 = tidyr::gather(post[c(T, rep(F, nSkip)), ], "anc", "p", 3:5) %>%
    filter(., p == max_p) %>% # only keep highest posterior prob. ancestry
    arrange(., chrom, position)
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
  post$max_p = apply(post[ , 3:5], 1, max)
  seg1 = tidyr::gather(post[c(T, rep(F, nSkip)), ], "anc", "p", 3:5) %>%
    filter(., p == max_p) %>% # only keep highest posterior prob. ancestry
    arrange(., chrom, position)
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
#maize = read.table(paste0(dir_in, "/maize_posteriors_small.txt"),
#                   sep = "\t", header = T, stringsAsFactors = F)
mex = do.call(rbind,
              lapply(include_mex$n, function(i) getPostSmall(id = i)))
#write.table(mex, paste0(dir_in, "/mex_posteriors_small.txt"),
#            sep = "\t", col.names = T, row.names = F,
#            quote = F)
#mex = read.table(paste0(dir_in, "/mex_posteriors_small.txt"),
#                   sep = "\t", header = T, stringsAsFactors = F)


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

# try to figure out the boundaries:
peak2_mex = do.call(rbind, 
                      lapply(include_mex$n, 
                             function(i) 
                               getPostSeg(id = i, chr = 4, start = 23000000,
                                          end = 28000000, buffer = 1000000,
                                          nSkip = 9))) # every 10th point plot

peak2_maize_plot <- plot_posteriors(post = peak2_maize, chr = 4, 
                start = 0, end = 1000000000,
                title = "maize_chr_4_peak2_zoom", save = F)
peak2_mex_plot <- plot_posteriors(post = peak2_mex, chr = 4, 
                                    start = 0, end = 1000000000,
                                    title = "mex_chr_4_peak2_zoom", save = F)

# I'll define peak2 on chr4 using this narrower region
#4:25800000-26400000
plot(peak2_maize_plot + 
       geom_vline(aes(xintercept=25900000)) +
       geom_vline(aes(xintercept=26300000)))
plot(peak2_mex_plot + 
       geom_vline(aes(xintercept=25800000)) +
       geom_vline(aes(xintercept=26400000)))
# show ZAnc peak?
# Now I classify individuals as maize_0, mex_1 etc. 
# where 0,1,2 is # mexicana alleles
head(peak2_mex)
peak2_anc <- rbind(peak2_mex, peak2_maize) %>% 
  filter(chrom == 4 & position > 25800000 & position < 26400000) %>%
  count(hilo_id, anc) %>%
  group_by(hilo_id) %>%
  slice(which.max(n)) %>%
  rename(top_anc_n = n) %>%
  left_join(., pass1[ , c("ID", "n", "zea")], by = c("hilo_id" = "n"))

table(peak2_anc$top_anc_n)
table(peak2_anc[,c("zea", "anc", "top_anc_n")])
# say 20/25 or more of one category can be classified
peak2_anc$anc[peak2_anc$top_anc_n<20] <- "unc.unc" # unknown

# write files to folder for outlier region
peak2_name <- "chr4/peak2"
peak2_dir <- paste0("../data/outliers/", 
                    peak2_name)
dir.create(peak2_dir,
           recursive = T)
peak2_anc <- peak2_anc %>% 
  mutate(hap_group = paste(zea, anc, sep = "_"))
peak2_hap_groups <- c(unique(peak2_anc$hap_group), allo_haps)
peak2_hap_files <- c(paste0("outliers/", peak2_name,"/", unique(peak2_anc$hap_group), ".list"), 
                     paste0("pass1_bam_pops/", allo_haps, ".list"))
# create all pairs of haplotype groups
peak2_hap_groups_pairs <- t(combn(peak2_hap_groups, m=2))
#peak2_hap_files_pairs <- t(combn(peak2_hap_files, m=2))
# write pairs to a file
for (i in 1:2){
  write.table(x = peak2_hap_groups_pairs[ , i], 
              file = paste0(peak2_dir, "/", "hap_groups_pair", i, ".list"),
              col.names = F, row.names = F, quote = F) 
}

# write a file with all haplotype groups that exist in this focal locus
write.table(x = peak2_hap_groups, 
            file = paste0(peak2_dir, "/", "hap_groups.list"),
            row.names = F, col.names = F, quote = F) 
# write file list
write.table(x = peak2_hap_files, 
            file = paste0(peak2_dir, "/", "hap_groups.files"),
            col.names = F, row.names = F, quote = F) 

# write a file linking each hilo ID to its hapgroup:
write.table(x = peak2_anc,
            file = paste0(peak2_dir, "/", "hap_groups.txt"),
            col.names = T, row.name = F, quote = F, sep = "\t")


# write a bam file list for each haplotype group
for (i in unique(peak2_anc$hap_group)){
  write.table(x = sapply(peak2_anc[peak2_anc$hap_group == i, "hilo_id"],
                         function(x) paste0("filtered_bam/hilo_", x, 
                         ".sort.dedup.baq.bam")),
              file = paste0(peak2_dir, "/", i, ".list"),
              col.names = F, row.names = F, quote = F)
}
# write regions file
write.table(x = "4:25800000-26400000",
            file = paste0(peak2_dir, "/", "regions.txt"),
            col.names = F, row.names = F, quote = F)

# a different peak - inversion on chr4: Inv4m
inv4m_post = do.call(rbind, 
                    lapply(c(include_mex$n, include_maize$n), 
                           function(i) 
                             getPostSeg(id = i, chr = 4, 
                                        start = inv[inv$ID=="inv4m" & inv$chr == 4, "start"],
                                        end = inv[inv$ID=="inv4m" & inv$chr == 4, "end"], 
                                        buffer = 0,
                                        nSkip = 99))) # don't need all the points
inv4m_anc <- inv4m_post %>% 
  count(hilo_id, anc) %>%
  group_by(hilo_id) %>%
  slice(which.max(n)) %>%
  rename(top_anc_n = n) %>%
  left_join(., pass1[ , c("ID", "n", "zea")], by = c("hilo_id" = "n"))

table(inv4m_anc$top_anc_n)
table(inv4m_anc[,c("zea", "anc", "top_anc_n")])
# say 80% or more of one category can be classified
inv4m_anc$anc[inv4m_anc$top_anc_n<.8*max(inv4m_anc$top_anc_n)] <- "unc.unc" # unknown

# write files to folder for outlier region
inv4m_name <- "chr4/inv4m"
inv4m_dir <- paste0("../data/outliers/", 
                    inv4m_name)
dir.create(inv4m_dir,
           recursive = T)
inv4m_anc <- inv4m_anc %>% 
  mutate(hap_group = paste(zea, anc, sep = "_"))
# write a file with all haplotype groups that exist in this focal locus
inv4m_hap_groups <- c(unique(inv4m_anc$hap_group), allo_haps)
inv4m_hap_files <- c(paste0("outliers/", inv4m_name,"/", unique(inv4m_anc$hap_group), ".list"), 
                     paste0("pass1_bam_pops/", allo_haps, ".list"))
# create all pairs of haplotype groups
inv4m_hap_groups_pairs <- t(combn(inv4m_hap_groups, m=2))
#inv4m_hap_files_pairs <- t(combn(inv4m_hap_files, m=2))
# write pairs to a file
for (i in 1:2){
  write.table(x = inv4m_hap_groups_pairs[ , i], 
              file = paste0(inv4m_dir, "/", "hap_groups_pair", i, ".list"),
              col.names = F, row.names = F, quote = F) 
}

# write all haplotype groups to a file
write.table(x = inv4m_hap_groups, 
            file = paste0(inv4m_dir, "/", "hap_groups.list"),
            col.names = F, row.names = F, quote = F) 
# write a file with all haplotype groups files
write.table(x = inv4m_hap_files, 
            file = paste0(inv4m_dir, "/", "hap_groups.files"),
            col.names = F, row.names = F, quote = F) 
# write a file linking each hilo ID to its hapgroup:
write.table(x = inv4m_anc,
            file = paste0(inv4m_dir, "/", "hap_groups.txt"),
            col.names = T, row.name = F, quote = F, sep = "\t")

# write a bam file for each
for (i in unique(inv4m_anc$hap_group)){
  write.table(x = sapply(inv4m_anc[inv4m_anc$hap_group == i, "hilo_id"],
                         function(x) paste0("filtered_bam/hilo_", x, 
                                            ".sort.dedup.baq.bam")),
              file = paste0(inv4m_dir, "/", i, ".list"),
              col.names = F, row.names = F, quote = F)
}
# write regions file
write.table(x = "4:171771502-185951149",
            file = paste0(inv4m_dir, "/", "regions.txt"),
            col.names = F, row.names = F, quote = F)

# now with these files I can calculate Fst between all the haplotype groups
# and pi within, and even plot the local ancestry tracts
# in PCA







# plot all of the inversions
inv = read.table("../data/refMaize/inversions/knownInv_v4_coord.txt",
                 stringsAsFactors = F, header = F)
colnames(inv) = c("ID", "chr", "start", "end", "length")

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

plot_inversion(post = rbind(maize, mex), 
               chr = inv$chrom[1], 
               start = inv$start[1], 
               end = inv$end[1], 
               ID = paste(inv$ID[1]), save = F)



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



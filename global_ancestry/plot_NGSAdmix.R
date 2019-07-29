# calculation/plotting global ancestry for low-coverage maize/mexicana individuals 
# in HILO adaptation project
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)# for color palette
# standard colors and labels for plots
blues = brewer.pal(n = 9, name = "Blues")[c(6,8)]
yellows = brewer.pal(n = 9, name = "YlOrBr")[c(4,6)]
colors_maize2mex = c(yellows, blues)
labels_maize2mex = c("allopatric_maize", "sympatric_maize", "sympatric_mexicana", "allopatric_mexicana")
colorsK = c(colors_maize2mex[c(4,1)], brewer.pal(n = 9, name = "YlOrRd")[6], brewer.pal(n=9, name = "YlGn")[6])

# minimum coverage to include in ancestry_hmm
min_coverage = 0.05

K = 2 # K = number of genetic clusters/groups
K = 3
K = 4

ancestries <- c("maize", "mexicana", "parviglumis")[1:K]

# starting with pass1 analysis from 1st round of sequencing
#file_prefix = paste0("../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedPCA/NGSAdmix/K", K)
PREFIX_METRICS <- "hilo_alloMAIZE_MAIZE4LOW"
PREFIX <- "pass2_alloMAIZE"
file_prefix <- paste0("results/NGSAdmix/", PREFIX, "/K", K)
admix <- read.table(paste0(file_prefix, ".qopt"))
colnames(admix) <- paste0("anc", 1:K) #c("anc1", "anc2")
out_dir <- paste0("results/NGSAdmix/", PREFIX)

# get metadata for individuals included in NGSadmix analysis
#pass1_allo4Low <- read.table("../data/pass1_allo4Low_ids.txt", stringsAsFactors = F, 
#                             header = T, sep = "\t")

IDs <- data.frame(ID = read.table(paste0("../samples/", PREFIX, "_IDs.list"), header = F, stringsAsFactors = F)$V1, stringsAsFactors = F)
metrics <- read.table(paste0("../filtered_bams/metrics/", PREFIX_METRICS, ".flagstat.total"), header = F, stringsAsFactors = F)
colnames(metrics) <- c("ID", "total_reads_pass")
hilo <- read.table("../samples/hilo_meta.txt", stringsAsFactors = F, header = T, sep = "\t")
maize4low <- read.table("../samples/MAIZE4LOW_meta.txt", stringsAsFactors = F, header = T)
landraces <- read.table("../samples/alloMAIZE_meta.txt", stringsAsFactors = F, header = T)
seq_Jan2019 <- read.table("../data/HILO_raw_reads/Jan2019_IDs.list", stringsAsFactors = F, header = F)$V1
parviglumis <- read.table("../samples/PalmarChico_meta.txt", stringsAsFactors = F, header = T, sep = "\t")


# only for using RIMMA0625_small
#landraces$ID[landraces$ID == "RIMMA0625"] <- "RIMMA0625_small"

# a rough coverage estimate is number of reads passing filtering * 150bp/read
# divided by total area they could map to in maize reference genome v4

# size of reference genome reads are mapped to
ref_genome_size <- sum(read.table("../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.fai", 
                                  stringsAsFactors = F)$V2) # V2 is the size of each chromosome mapped to,
# including parts I won't analyze on the Pt and Mt and scaffolds not assigned to chromosomes (but reads would pass Q filters there too)
metrics$est_coverage = round(metrics$total_reads_pass*150/ref_genome_size, 4)

# get elevation data
elev = read.table("../data/riplasm/gps_and_elevation_for_sample_sites.txt",
                  stringsAsFactors = F, header = T, sep = "\t")


# combine sample meta data
meta <- bind_rows(hilo, maize4low, landraces, parviglumis) %>%
  left_join(., metrics, by = "ID") %>%
  mutate(., group = paste(symp_allo, zea, sep = "_")) %>%
  left_join(., elev, by = c("popN", "zea", "symp_allo", 
                            "RI_ACCESSION", "GEOCTY", "LOCALITY")) %>%
  mutate(est_coverage = ifelse(zea == "parviglumis", 10, est_coverage))

# join bams and admix by position (CAUTION - bam list order and admix results MUST MATCH!)
d <- bind_cols(IDs, admix)  %>%
  left_join(., meta, by = "ID") %>%
  arrange(., popN) %>%
  arrange(., zea) %>%
  arrange(., symp_allo)

# order by elevation and group, filter out low coverage ind's
d_meta <- d %>%
  filter(., est_coverage >= min_coverage) %>% # only include individuals with at least .05x coverage
  .[with(., order(group, ELEVATION)), ]

# for K=2 only, which ancestry anc1 or anc2 is mexicana-like?
if (K == 2) {
  mex_anc  <- d %>%
    filter(group == "allopatric_mexicana") %>%
    dplyr::select(c("anc1", "anc2")) %>%
    apply(., 2, mean) %>%
    which.max(.)
  maize_anc  <- d %>%
    filter(group == "allopatric_maize") %>%
    dplyr::select(c("anc1", "anc2")) %>%
    apply(., 2, mean) %>%
    which.max(.)
  d2 <- dplyr::mutate(d, alpha_mex = round(d[ , names(mex_anc)], 10)) %>%
    dplyr::mutate(., alpha_maize = round(d[ , names(maize_anc)], 10))

  # make anc_hmm input file with population #, avg. maize proportion, avg. mex proportion ancestry
  # first make directory
  #path_global_out = file.path("../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm", "input")
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = T)
  # then create data summarizing mean global ancestry for each pop
  # excluding ind's w/ low coverage. for Ancestry_hmm input
  alphasByPop = d2 %>%
    filter(., symp_allo == "sympatric" & est_coverage >= min_coverage) %>%
    group_by(., popN) %>%
    summarise(., alpha_maize = mean(alpha_maize), 
              alpha_mex = mean(alpha_mex),
              n = n()) # number included individuals per group
  table(filter(d, symp_allo == "sympatric" & est_coverage >= min_coverage)$popN) # numbers of individuals per pop
  #View(cbind(alphasByPop,table(filter(d, symp_allo == "sympatric" & est_coverage >= 0.05)$popN) ))
  # write file with mean population ancestry (for prior in ancestry_hmm)
  write.table(format(as.data.frame(alphasByPop), scientific = F),  # write table so I can pull population global ancestry values for ancestry_hmm
            file.path(out_dir, "globalAdmixtureIncludedByPopN.txt"),
            col.names = T, row.names = F, quote = F, sep = "\t")
  # also write out individual admixture proportions for included individuals
  ind_included_summary <- d2 %>%
    filter(., symp_allo == "sympatric" & est_coverage >= min_coverage) %>%
    arrange(., n) %>%
    arrange(., popN) %>%
    dplyr::select(., c("ID", "popN", "alpha_mex", "alpha_maize", "est_coverage"))
  
  write.table(x = format(ind_included_summary, digits = 10, scientific = F),  # write table so I can pull population global ancestry values for ancestry_hmm
              file = file.path(out_dir, "globalAdmixtureByIncludedIndividual.txt"),
              col.names = T, row.names = F, quote = F, sep = "\t")
}

# lets look at allele frequency estimates for each ancestry
#freqs <- read.table("../data/NGSadmix/pass1/K2_pruned_all.fopt.gz") 
freqs <- read.table(paste0(file_prefix, ".fopt.gz"))
colnames(freqs) <- paste0("anc", 1:K)
#freqs$diff <- abs(freqs$anc1 - freqs$anc2)
mean_allele_freq <- apply(freqs, 1, mean)
hetSubpops <- sapply(1:K, function(i) 
  mean(2 * freqs[ , paste0("anc", i)] * 
         (1 - freqs[ , paste0("anc", i)])))
#het1 <- mean(freqs$anc1 * (1-freqs$anc1))
#het2 <- mean(freqs$anc2 * (1-freqs$anc2))
# a better total approximate would be to use the MAF of the entire sample estimated in ANGSD
# rather than the mean of equal sized hypothetical pop 'ancestries'
hetTot <- mean(2 * mean_allele_freq * (1 - mean_allele_freq))
# differentiation is fairly low between the two identified genetic clusters
Fst <- (hetTot - mean(hetSubpops))/hetTot
Fst

# make STRUCTURE-like ancestry plots:
p1 <- d %>% 
  filter(est_coverage > min_coverage) %>%
  arrange(., popN) %>%
  tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  ggplot(., aes(fill=ancestry, y=p, x=ID)) +
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~group) +
  ggtitle(paste0("pass 2 K=2 NGSAdmix results. ind's with depth > ", min_coverage))
plot(p1)
ggsave(paste0("plots/NGS_admix_K", K, "_", PREFIX, ".png"), 
       plot = p1, 
       device = "png", 
       width = 15, height = 8, units = "in",
       dpi = 200)
# plot just allopatric maize
p1_allomaize <- d %>% 
  filter(est_coverage > 0.05) %>%
  filter(group == "allopatric_maize") %>%
  tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  ggplot(., aes(fill=ancestry, y=p, x=ID)) +
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~LOCALITY) +
  ggtitle(paste0("pass 2 K=2 NGSAdmix results. ind's with depth > ", min_coverage))
plot(p1_allomaize)
ggsave(paste0("plots/NGS_admix_allo_maize_only_K", K, "_", PREFIX, ".png"), 
       plot = p1_allomaize, 
       device = "png", 
       width = 15, height = 8, units = "in",
       dpi = 200)



# plot 'STRUCTURE-like' ancestry plots
png(paste0("../plots/NGSadmix_K", K, ".png"), # saves plot as pin in ../plots/
    height = 5, width = 8, units = "in", res = 150)
d %>%
  dplyr::select(., colnames(admix)) %>%
  t(.) %>%
  barplot(height = .,
          col = colorsK[1:K],
          space=0,
          border=NA,
          main = "all HILO populations",
          xlab="Individuals",
          ylab="admixture") 
mtext(paste("est. between-ancestry Fst:", round(Fst,2)))
dev.off()

# now plot again with individuals with very low coverage <0.05x filtered out      
png(paste0("../plots/NGSadmix_K", K, "_over_", min_coverage, "x_coverage.png"),
    height = 5, width = 8, units = "in", res = 150)
d %>%
  filter(., est_coverage >= min_coverage) %>%
  select(., colnames(admix)) %>%
  t(.) %>%
  barplot(height = .,
          col = colorsK[1:K],
          space=0,
          border=NA,
          main = "all HILO populations",
          xlab=paste0("Individuals with est. coverage >= ", min_coverage, "x"),
          ylab="admixture") 
mtext(paste("est. between-ancestry Fst:", round(Fst,2)))
dev.off()

for (g in unique(d$group)){
  png(paste0("../plots/NGSadmix_K", K, "_", g, ".png"), # saves plot as pin in ../plots/
      height = 5, width = 8, units = "in", res = 150)
  d %>%
    filter(., group == g) %>%
    filter(., est_coverage >= min_coverage) %>%
    select(., colnames(admix)) %>%
    t(.) %>%
    barplot(height = .,
            col = colorsK[1:K],
            space=0,
            border=NA,
            main = g,
            xlab=paste0("Individuals with est. coverage >= ", min_coverage, "x"),
            ylab="admixture")  
  dev.off()
}


# by location
if (K == 2){
  p = d %>%
    filter(., est_coverage >= .05) %>%
    mutate(group = factor(paste(symp_allo, zea, sep = "_"), ordered = T,
                          levels = labels_maize2mex)) %>%
    ggplot(., aes(x = group, y = anc1, color = group)) +
    geom_violin() + 
    geom_point(size=1, position = position_jitter(w=0.05)) +
    scale_colour_manual(values = colors_maize2mex) +
    #scale_colour_manual(values = c("yellow", "darkblue", "orange", "blue")) +
    theme_minimal() + 
    labs(y = "proportion 'mexicana-like' ancestry")
  p + facet_wrap(~LOCALITY) # show plot
  # save plot
  ggsave("plots/NGSadmix_proportion_mexicana-like_by_location.png", plot = p + facet_wrap(~LOCALITY), device = png(), 
         width = 12, height = 8, units = "in",
         dpi = 200)
}


# make NGSadmix plot for poster:
# get elevation data
meta.pops = read.table("../data/riplasm/gps_and_elevation_for_sample_sites.txt",
                 stringsAsFactors = F, header = T, sep = "\t")
meta.ind = left_join(d[ , c("ID", "popN", "total_reads_pass", "est_coverage", "group")], meta.pops, by = "popN") %>%
  filter(., est_coverage >= min_coverage) %>% # only include individuals with at least .05x coverage
  .[with(., order(group, ELEVATION)), ]
# combined with admixture results
d_meta.ind = left_join(dplyr::select(d, ID, starts_with("anc")), meta.ind, by = "ID")


# simple linear model predicting ancestry2 from group and elevation
d_meta.ind %>%
  filter(., symp_allo == "sympatric") %>%
  lm(data = ., anc2 ~ zea + ELEVATION + est_coverage) %>%
  summary(.)
lmZeaElev <- d_meta.ind %>%
  filter(., symp_allo == "sympatric") %>%
  lm(data = ., anc2 ~ zea*ELEVATION) 
summary(lmZeaElev)
# plot linear model predictions for effects of elevation
d_meta.ind %>%
  filter(., symp_allo == "sympatric") %>%
  ggplot(., aes(x = ELEVATION, y = anc2, color = zea, size = est_coverage)) +
  geom_point() +
  ylab("mexicana ancestry") +
  scale_colour_manual(values = colors_maize2mex[2:3]) +
  geom_abline(intercept = lmZeaElev$coefficients["(Intercept)"] + lmZeaElev$coefficients["zeamexicana"], 
              slope = lmZeaElev$coefficients["ELEVATION"] + lmZeaElev$coefficients["zeamexicana:ELEVATION"],
              color = colors_maize2mex[3]) +
  geom_abline(intercept = lmZeaElev$coefficients["(Intercept)"], 
              slope = lmZeaElev$coefficients["ELEVATION"],
              color = colors_maize2mex[2]) +
  ggtitle("Higher mexicana ancestry at higher elevations")
ggsave("plots/lm_predict_NGSadmix_proportion_mexicana-like_by_elevation.png", 
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)
d_meta.ind %>%
  filter(., symp_allo == "sympatric") %>%
  ggplot(., aes(x = ELEVATION, y = anc2, color = LOCALITY, size = est_coverage)) +
  geom_point() +
  ylab("mexicana ancestry") +
  #scale_colour_manual(values = colors_maize2mex[2:3]) +
  geom_abline(intercept = lmZeaElev$coefficients["(Intercept)"] + lmZeaElev$coefficients["zeamexicana"], 
              slope = lmZeaElev$coefficients["ELEVATION"] + lmZeaElev$coefficients["zeamexicana:ELEVATION"],
              color = colors_maize2mex[3]) +
  geom_abline(intercept = lmZeaElev$coefficients["(Intercept)"], 
              slope = lmZeaElev$coefficients["ELEVATION"],
              color = colors_maize2mex[2]) +
  ggtitle("Higher mexicana ancestry at higher elevations")
ggsave("plots/lm_predict_NGSadmix_proportion_mexicana-like_by_elevation_colored_by_pop.png", 
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)

# now plot again with individuals with very low coverage <0.05x filtered out      
png(paste0("plots/", PREFIX, "NGSadmix_K", K, "_over_", min_coverage, "x_coverage_wElevation_for_poster.png"),
    height = 5, width = 8, units = "in", res = 300)
par(mar=c(4.1,4.1,4.1,4.1))
bar = d_meta.ind %>%
  #dplyr::select(., colnames(admix)) %>%
  dplyr::select(., rev(colnames(admix))) %>% # reverse order so blue = mexicana
  t(.) %>%
  barplot(height = .,
          col = colorsK[1:K],
          space=0,
          border=NA, xaxt = "n")
title(main = "Admixture in highland maize and mexicana",
      ylab=paste0("Ancestry proportion K=", K, " (NGSAdmix)"))
#title(xlab="Allopatric Maize | Allopatric Mexicana | Sympatric Maize | Sympatric Mexicana", 
#      line = 1)
text(x = tapply(bar, d_meta.ind$group, mean), par("usr")[3], srt = 60, adj= 1, xpd = TRUE,
     labels = unique(d_meta.ind$group), cex=0.45)
#title(xlab=paste("Allopatric Ref.", "Sympatric Maize", "Sympatric Mexicana",
#      sep = "               |               "),
#      line = 1)
par(new = TRUE)
plot(x = bar, y = d_meta.ind$ELEVATION, ylim = c(1000, 3000), axes = F, type = "p", cex = .1, 
     xlab = "", ylab = "")
axis(side = 4, at = c(1000, 2000, 3000), 
     tick = T, labels = T, col = "black")
mtext("Elevation (m)", side=4, line=2.5)
par(mar=c(5.1,4.1,4.1,2.1)) # set back default
dev.off()

# make simple plots for SNPs pruned every 1000th SNP
file_prefix_1000 = paste0("../data/geno_lik/merged_pass1_all_alloMaize4Low_16/prunedBy1000/NGSAdmix/K", K)

admix_1000 <- read.table(paste0(file_prefix_1000, ".qopt"))
colnames(admix_1000) <- paste0("anc", 1:K) #c("anc1", "anc2")
admix_1000 <- admix_1000[, K:1] # reverse order for pruned by 1000
d_1000 <- bind_cols(pass1_allo4Low, admix_1000)  %>%
  arrange(., popN) %>%
  arrange(., zea) %>%
  arrange(., symp_allo) %>%
  mutate(., group = paste(symp_allo, zea, sep = "_"))
d_meta.ind_1000 = left_join(d_1000, meta, by = c("popN", "zea", "symp_allo", 
                                   "RI_ACCESSION", "GEOCTY", "LOCALITY")) %>%
  filter(., est_coverage >= .05) %>% # only include individuals with at least .05x coverage
  .[with(., order(group, ELEVATION)), ]


# plot results for SNPs pruned down to 1 per 1000 SNPs
png(paste0("../plots/NGSadmix_K", K, "_over_0.05x_coverage_wElevation_for_poster_pruned1000.png"),
    height = 5, width = 8, units = "in", res = 300)
par(mar=c(4.1,4.1,4.1,4.1))
bar_1000 = d_meta.ind_1000 %>%
  select(., colnames(admix_1000)) %>%
  t(.) %>%
  barplot(height = .,
          col = colorsK[1:K],
          space=0,
          border=NA, xaxt = "n")
title(main = "Admixture in highland maize and mexicana; every 1000th SNP",
      ylab=paste0("Ancestry proportion K=", K, " (NGSAdmix)"))
title(xlab="Allopatric Maize | Allopatric Mexicana | Sympatric Maize | Sympatric Mexicana", 
      line = 1)

par(new = TRUE)
plot(x = bar_1000, y = d_meta.ind$ELEVATION, ylim = c(1000, 3000), axes = F, type = "p", cex = .1, 
     xlab = "", ylab = "")
axis(side = 4, at = c(1000, 2000, 3000), 
     tick = T, labels = T, col = "black")
mtext("Elevation (m)", side=4, line=2.5)
par(mar=c(5.1,4.1,4.1,2.1)) # set back default
dev.off()

# do global ancestry estimatesiloveyou
# from NGSadmix confirm the ancestry ~ r pattern?
# get ngsadmix k=2 estimates for 5 bins of recombination rate:
admix_r <- do.call(rbind,
                   lapply(1:5, function(i)
                     read.table(paste0("results/NGSAdmix/", PREFIX, "/recomb_", i,"_1percent/K2.qopt")) %>%
                       bind_cols(IDs, .) %>%
                       left_join(., meta.ind, by = "ID") %>%
                       arrange(., popN) %>%
                       arrange(., zea) %>%
                       arrange(., symp_allo) %>%
                       filter(est_coverage > .05) %>%
                       dplyr::mutate(., recomb_bin = paste0('recomb_', i)) %>%
                       dplyr::mutate(., # label 'mexicana' ancestry least common in allopatric maize
                                     mexicana_ancestry = sapply(1:nrow(.), function(j) ifelse(mean(filter(., group == "allopatric_maize")$V1) < .5, 
                                                                                              .$V1[j], 
                                                                                              .$V2[j])))
                   ))
# plot across recombination bins:
admix_r %>%
  dplyr::mutate(est_coverage = ifelse(est_coverage > 2, 2, est_coverage)) %>%
  filter(est_coverage > .05) %>%
  ggplot(aes(x = recomb_bin, y = mexicana_ancestry, 
             color = LOCALITY, size = est_coverage)) +
  geom_point() +
  facet_wrap(~group) +
  ggtitle("individual K=2 NGSAdmix mex-like ancestry estimate by r bin") +
  xlab("low to high recombination bins (quintiles of 10kb windows)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/ind_mexicana-like-ancestry_by_recomb_bin.png",
       height = 8, width = 8,
       units = "in", device = "png")

# plot across recombination bins:
admix_r %>%
  dplyr::mutate(est_coverage = ifelse(est_coverage > 2, 2, est_coverage)) %>%
  filter(est_coverage > 0.05) %>%
  ggplot(aes(x = recomb_bin, y = mexicana_ancestry, 
             color = LOCALITY, group = ID)) +
  geom_line() +
  #geom_point(aes(size = est_coverage)) +
  facet_wrap(~group) +
  ggtitle("individual K=2 NGSAdmix mex-like ancestry estimate by r bin") +
  xlab("low to high recombination bins (quintiles of 10kb windows)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/ind_mexicana-like-ancestry_by_recomb_bin_lines.png",
       height = 8, width = 8,
       units = "in", device = "png")

# does the strength of selection against mexicana ancestry seem stronger 
# for lower elevation populations?
admix_r %>%
  dplyr::mutate(est_coverage = ifelse(est_coverage > 2, 2, est_coverage)) %>%
  filter(est_coverage > 0.05) %>%
  ggplot(aes(x = recomb_bin, y = mexicana_ancestry, 
             color = ELEVATION, group = ID)) +
  geom_line() +
  #geom_point(aes(size = est_coverage)) +
  facet_wrap(~group) +
  ggtitle("individual K=2 NGSAdmix mex-like ancestry estimate by r bin") +
  xlab("low to high recombination bins (quintiles of 10kb windows)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/ind_mexicana-like-ancestry_by_recomb_bin_lines_colorByElevation.png",
       height = 8, width = 8,
       units = "in", device = "png")

# separate pattern by population -- seems to hold across most populations individually
admix_r %>%
  filter(., est_coverage > 0.05) %>%
  ggplot(aes(fill = recomb_bin, x = reorder(LOCALITY, ELEVATION), y = mexicana_ancestry)) +
  geom_boxplot() +
  facet_wrap(~group) +
  ggtitle("NGSadmix population mex-like ancestry estimates by r bin (boxplot of ind admixture proportions)") +
  xlab("low to high recombination bins (quintiles of 10kb windows)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/pop_mexicana-like-ancestry_by_recomb_bin_pops_by_locality.png",
       height = 10, width = 12,
       units = "in", device = "png")

admix_r %>%
  filter(., est_coverage > 0.05) %>%
  ggplot(aes(fill = recomb_bin, y = mexicana_ancestry)) +
  geom_boxplot() +
  facet_wrap(~group) +
  ggtitle("NGSadmix mex-like ancestry estimates by r bin (boxplot of ind admixture proportions)") +
  xlab("low to high recombination bins (quintiles of 10kb windows)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/pop_mexicana-like-ancestry_by_recomb_bin_boxplot.png",
       height = 6, width = 8,
       units = "in", device = "png")
# hmm .. looks like everyone has more mexicana ancestry in high r windows
# including mexicana

# remake plot of ancestry ~ elevation for different recombination bins
admix_r %>%
  filter(., symp_allo == "sympatric") %>%
  ggplot(., aes(x = ELEVATION, y = mexicana_ancestry, color = zea, size = est_coverage)) +
  geom_point() +
  ylab("mexicana ancestry") +
  scale_colour_manual(values = colors_maize2mex[2:3]) +
  geom_smooth(method = "lm", aes(group = zea)) +
  ggtitle("Higher mexicana ancestry at higher elevations") +
  facet_wrap(~recomb_bin)
ggsave("plots/lm_predict_NGSadmix_proportion_mexicana-like_by_elevation_andByRecombBin.png", 
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)

# color this plot by LOCALITY
admix_r %>%
  filter(., symp_allo == "sympatric") %>%
  ggplot(., aes(x = ELEVATION, y = mexicana_ancestry, color = LOCALITY, size = est_coverage)) +
  geom_point() +
  ylab("mexicana ancestry") +
  geom_smooth(method = "lm", aes(group = zea), color = "black") +
  ggtitle("Higher mexicana ancestry at higher elevations") +
  facet_wrap(~recomb_bin)
ggsave("plots/lm_predict_NGSadmix_proportion_mexicana-like_by_elevation_andByRecombBin_colorByPop.png", 
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)





# now look at Fst for the different r bins
freqs_f <- do.call(rbind,
                   lapply(1:5, function(i)
                     read.table(paste0("results/NGSAdmix/", PREFIX, "/recomb_", i,"_1percent/K2.fopt.gz")) %>%
                       dplyr::mutate(., het1 = 2*V1*(1-V1)) %>%
                       dplyr::mutate(., het2 = 2*V2*(1-V2)) %>%
                       dplyr::mutate(., pTot = (V1+V2)/2) %>%
                       dplyr::mutate(., hetTot = 2*pTot*(1-pTot)) %>%
                       dplyr::mutate(., piBetween = V1*(1-V2) + (1-V1)*V2) %>%
                       dplyr::mutate(., recomb_bin = paste0("recomb_", i))
                   ))

maize_is_anc1 <- admix_r %>%
  filter(., group == "allopatric_maize") %>%
  group_by(recomb_bin) %>%
  dplyr::summarise(isTrue = mean(V1) > .5)
                                                                         
freqs_f %>%
  dplyr::group_by(recomb_bin) %>%
  dplyr::summarise(., Fst = 1 - (mean(het1) + mean(het2))/2/mean(hetTot),
                   het1 = mean(het1),
                   het2 = mean(het2),
                   hetTot = mean(hetTot),
                   piBetween = mean(piBetween)) %>%
  left_join(., maize_is_anc1, by = "recomb_bin") %>%
  dplyr::mutate(hetMaizeAnc = ifelse(maize_is_anc1$isTrue, het1, het2)) %>%
  dplyr::mutate(hetMexicanaAnc = ifelse(maize_is_anc1$isTrue, het2, het1)) %>%
  gather(., "stat", "value", c("hetMaizeAnc", "hetMexicanaAnc", "piBetween", "Fst", "hetTot")) %>%
  ggplot(aes(x = recomb_bin, y = value, group = stat, color = stat)) +
  geom_point() +
  geom_line() +
  facet_wrap(~stat=="Fst", scales = "free") +
  ggtitle("pi and Fst for ancestry-estimated-allele-freqs from NGSadmix K=2 clusters") +
  xlab("low to high recomb bins (quintiles of 10kb windows)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/fst_by_recomb_bin.png",
       height = 6, width = 12,
       units = "in", device = "png")

# make plots for pruned every 100th SNP:
K = 2
PREFIX = "pass2_alloMAIZE"
IDs <- data.frame(ID = read.table(paste0("../samples/", PREFIX, "_IDs.list"), header = F, stringsAsFactors = F)$V1, stringsAsFactors = F)
ancestries <- c("maize", "mexicana", "parviglumis")[1:K]


file_prefix_by100 <- paste0("results/NGSAdmix/", PREFIX, "/prunedBy100/K", K)
admix_by100 <- read.table(paste0(file_prefix_by100, ".qopt"))
colnames(admix_by100) <- paste0("anc", 1:K) #c("anc1", "anc2")
d_by100 <- bind_cols(IDs, admix_by100)  %>%
  left_join(., meta, by = "ID") %>%
  arrange(., popN) %>%
  arrange(., zea) %>%
  arrange(., symp_allo) %>%
  filter(., est_coverage >= .05) %>% # only include individuals with at least .05x coverage
  .[with(., order(group, ELEVATION)), ]

# assign mexicana ancestry
which_anc_by100 <- data.frame(ancestry = colnames(admix_by100),
                        ancestry_label = 
                          sapply(colnames(admix_by100), function(x) 
                            names(which.max(tapply(d_by100[ , x], d_by100$zea, mean)))),
                        stringsAsFactors = F)

anc_by100 = d_by100 %>% 
  tidyr::gather(., "ancestry", "p", colnames(admix_by100)) %>%
  left_join(., which_anc_by100, by = "ancestry") %>%
  dplyr::select(-ancestry) %>%
  tidyr::spread(., ancestry_label, p)

# simple admixture plot:
anc_by100 %>% 
  arrange(., popN) %>%
  tidyr::gather(., "ancestry", "p", ancestries) %>%
  ggplot(., aes(fill=ancestry, y=p, x=ID)) +
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~group)

# plot linear models of ancestry over elevation
lmZeaElev_by100 <- anc_by100 %>%
  filter(., symp_allo == "sympatric") %>%
  lm(data = ., mexicana ~ zea*ELEVATION)

anc_by100 %>%
  filter(., symp_allo == "sympatric") %>%
  ggplot(., aes(x = ELEVATION, y = mexicana, color = LOCALITY, size = est_coverage, shape = zea)) +
  geom_point() +
  ylab("mexicana ancestry") +
  geom_abline(intercept = lmZeaElev_by100$coefficients["(Intercept)"] + lmZeaElev_by100$coefficients["zeamexicana"], 
              slope = lmZeaElev_by100$coefficients["ELEVATION"] + lmZeaElev_by100$coefficients["zeamexicana:ELEVATION"],
              color = colors_maize2mex[3]) +
  geom_abline(intercept = lmZeaElev_by100$coefficients["(Intercept)"], 
              slope = lmZeaElev_by100$coefficients["ELEVATION"],
              color = colors_maize2mex[2]) +
  ggtitle("Clines in mexicana ancestry across elevation") +
  #theme(legend.position="bottom") +
  labs(color = "Location", shape = "Subspecies")
ggsave(paste0("plots/lm_predict_NGSadmix_proportion_mexicana-like_by_elevation_colored_by_pop_prunedBy100_", PREFIX, "_K", K, ".png"), 
       device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)




        
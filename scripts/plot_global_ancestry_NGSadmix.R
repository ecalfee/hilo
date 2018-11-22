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

K = 2 # K = number of genetic clusters/groups
K = 3
K = 4
# starting with pass1 analysis from 1st round of sequencing
file_prefix = paste0("../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedPCA/NGSAdmix/K", K)
admix <- read.table(paste0(file_prefix, ".qopt"))
colnames(admix) <- paste0("anc", 1:K) #c("anc1", "anc2")

# get metadata for individuals included in NGSadmix analysis
pass1_allo4Low <- read.table("../data/pass1_allo4Low_ids.txt", stringsAsFactors = F, 
                             header = T, sep = "\t")
# join bams and admix by position (CAUTION - bam list order and admix results MUST MATCH!)
d <- bind_cols(pass1_allo4Low, admix)  %>%
  arrange(., popN) %>%
  arrange(., zea) %>%
  arrange(., symp_allo) %>%
  mutate(., group = paste(symp_allo, zea, sep = "_"))

# make anc_hmm input file with population #, avg. maize proportion, avg. mex proportion ancestry
# first make directory
path_global_out = file.path("../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm", "input")
if (!dir.exists(path_global_out)) dir.create(path_global_out, recursive = T)
# then create data summarizing mean global ancestry for each pop
# excluding ind's w/ low coverage. for Ancestry_hmm input
alphasByPop = d %>%
  filter(., symp_allo == "sympatric" & est_coverage >= 0.05) %>%
  group_by(., popN) %>%
  summarise(., alpha_maize = round(mean(anc2),2)) %>%
  mutate(., alpha_mex = 1 - alpha_maize)
table(filter(d, symp_allo == "sympatric" & est_coverage >= 0.05)$popN) # numbers of individuals per pop
#View(cbind(alphasByPop,table(filter(d, symp_allo == "sympatric" & est_coverage >= 0.05)$popN) ))
if (K == 2) {
  write.table(format(as.data.frame(alphasByPop), digits=2, scientific = F),  # write table so I can pull population global ancestry values for ancestry_hmm
            file.path(path_global_out, "globalAdmixtureByPopN.txt"),
            col.names = F, row.names = F, quote = F, sep = "\t")
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

# plot 'STRUCTURE-like' ancestry plots
png(paste0("../plots/NGSadmix_K", K, ".png"), # saves plot as pin in ../plots/
    height = 5, width = 8, units = "in", res = 150)
d %>%
  select(., colnames(admix)) %>%
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
png(paste0("../plots/NGSadmix_K", K, "_over_0.05x_coverage.png"),
    height = 5, width = 8, units = "in", res = 150)
d %>%
  filter(., est_coverage >= .05) %>%
  select(., colnames(admix)) %>%
  t(.) %>%
  barplot(height = .,
          col = colorsK[1:K],
          space=0,
          border=NA,
          main = "all HILO populations",
          xlab="Individuals with est. coverage >= 0.05x",
          ylab="admixture") 
mtext(paste("est. between-ancestry Fst:", round(Fst,2)))
dev.off()

for (g in unique(d$group)){
  png(paste0("../plots/NGSadmix_K", K, "_", g, ".png"), # saves plot as pin in ../plots/
      height = 5, width = 8, units = "in", res = 150)
  d %>%
    filter(., group == g) %>%
    filter(., est_coverage >= .05) %>%
    select(., colnames(admix)) %>%
    t(.) %>%
    barplot(height = .,
            col = colorsK[1:K],
            space=0,
            border=NA,
            main = g,
            xlab="Individuals with est. coverage >= 0.05x",
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
  ggsave("../plots/NGSadmix_proportion_mexicana-like_by_location.png", plot = p + facet_wrap(~LOCALITY), device = png(), 
         width = 12, height = 8, units = "in",
         dpi = 200)
}


# make NGSadmix plot for poster:
# get elevation data
meta = read.table("../data/riplasm/gps_and_elevation_for_sample_sites.txt",
                 stringsAsFactors = F, header = T, sep = "\t")
d_meta = left_join(d, meta, by = c("popN", "zea", "symp_allo", 
                                   "RI_ACCESSION", "GEOCTY", "LOCALITY")) %>%
  filter(., est_coverage >= .05) %>% # only include individuals with at least .05x coverage
  .[with(., order(group, ELEVATION)), ]


  # now plot again with individuals with very low coverage <0.05x filtered out      
png(paste0("../plots/NGSadmix_K", K, "_over_0.05x_coverage_wElevation_for_poster.png"),
    height = 5, width = 8, units = "in", res = 300)
par(mar=c(4.1,4.1,4.1,4.1))
bar = d_meta %>%
  select(., colnames(admix)) %>%
  t(.) %>%
  barplot(height = .,
          col = colorsK[1:K],
          space=0,
          border=NA, xaxt = "n")
title(main = "Admixture in highland maize and mexicana",
      ylab=paste0("Ancestry proportion K=", K, " (NGSAdmix)"))
title(xlab="Allopatric Maize | Allopatric Mexicana | Sympatric Maize | Sympatric Mexicana", 
      line = 1)
#title(xlab=paste("Allopatric Ref.", "Sympatric Maize", "Sympatric Mexicana",
#      sep = "               |               "),
#      line = 1)
par(new = TRUE)
plot(x = bar, y = d_meta$ELEVATION, ylim = c(1000, 3000), axes = F, type = "p", cex = .1, 
     xlab = "", ylab = "")
axis(side = 4, at = c(1000, 2000, 3000), 
     tick = T, labels = T, col = "black")
mtext("Elevation (m)", side=4, line=2.5)
par(mar=c(5.1,4.1,4.1,2.1)) # set back default
dev.off()

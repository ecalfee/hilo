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
#file_prefix = paste0("../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedPCA/NGSAdmix/K", K)
PREFIX <- "hilo_alloMAIZE_MAIZE4LOW_RIMMA0625_small"
file_prefix <- paste0("results/NGSAdmix/", PREFIX, "/K", K)
admix <- read.table(paste0(file_prefix, ".qopt"))
colnames(admix) <- paste0("anc", 1:K) #c("anc1", "anc2")

# get metadata for individuals included in NGSadmix analysis
#pass1_allo4Low <- read.table("../data/pass1_allo4Low_ids.txt", stringsAsFactors = F, 
#                             header = T, sep = "\t")

IDs <- data.frame(ID = read.table(paste0("../samples/", PREFIX, "_IDs.list"), header = F, stringsAsFactors = F)$V1, stringsAsFactors = F)
metrics <- read.table(paste0("../filtered_bams/metrics/", PREFIX, ".flagstat.total"), header = F, stringsAsFactors = F)
colnames(metrics) <- c("ID", "total_reads_pass")
hilo <- read.table("../samples/hilo_meta.txt", stringsAsFactors = F, header = T, sep = "\t")
maize4low <- read.table("../samples/MAIZE4LOW_meta.txt", stringsAsFactors = F, header = T)
landraces <- read.table("../samples/alloMAIZE_meta.txt", stringsAsFactors = F, header = T)
seq_Jan2019 <- read.table("../data/HILO_raw_reads/Jan2019_IDs.list", stringsAsFactors = F, header = F)$V1

# only for using RIMMA0625_small
landraces$ID[landraces$ID == "RIMMA0625"] <- "RIMMA0625_small"

# a rough coverage estimate is number of reads passing filtering * 150bp/read
# divided by total area they could map to in maize reference genome v4

# size of reference genome reads are mapped to
ref_genome_size <- sum(read.table("../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.fai", 
                                  stringsAsFactors = F)$V2) # V2 is the size of each chromosome mapped to,
# including parts I won't analyze on the Pt and Mt and scaffolds not assigned to chromosomes (but reads would pass Q filters there too)
metrics$est_coverage = round(metrics$total_reads_pass*150/ref_genome_size, 4)

# combine sample meta data
meta <- bind_rows(hilo, maize4low, landraces) %>%
  left_join(., metrics, by = "ID") %>%
  mutate(., group = paste(symp_allo, zea, sep = "_"))

# join bams and admix by position (CAUTION - bam list order and admix results MUST MATCH!)
d <- bind_cols(IDs, admix)  %>%
  left_join(., meta, by = "ID") %>%
  arrange(., popN) %>%
  arrange(., zea) %>%
  arrange(., symp_allo)

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

# make STRUCTURE-like ancestry plots:
p1 <- d %>% 
  filter(est_coverage > 0.05) %>%
  tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  ggplot(., aes(fill=ancestry, y=p, x=ID)) +
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~group) +
  ggtitle("pass 2 K=2 NGSAdmix results. ind's with depth > 0.05")
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
  ggtitle("pass 2 K=2 NGSAdmix results. ind's with depth > 0.05")
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
d_meta_1000 = left_join(d_1000, meta, by = c("popN", "zea", "symp_allo", 
                                   "RI_ACCESSION", "GEOCTY", "LOCALITY")) %>%
  filter(., est_coverage >= .05) %>% # only include individuals with at least .05x coverage
  .[with(., order(group, ELEVATION)), ]


# plot results for SNPs pruned down to 1 per 1000 SNPs
png(paste0("../plots/NGSadmix_K", K, "_over_0.05x_coverage_wElevation_for_poster_pruned1000.png"),
    height = 5, width = 8, units = "in", res = 300)
par(mar=c(4.1,4.1,4.1,4.1))
bar_1000 = d_meta_1000 %>%
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
plot(x = bar_1000, y = d_meta$ELEVATION, ylim = c(1000, 3000), axes = F, type = "p", cex = .1, 
     xlab = "", ylab = "")
axis(side = 4, at = c(1000, 2000, 3000), 
     tick = T, labels = T, col = "black")
mtext("Elevation (m)", side=4, line=2.5)
par(mar=c(5.1,4.1,4.1,2.1)) # set back default
dev.off()

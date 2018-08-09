# calculation/plotting global ancestry for low-coverage maize/mexicana individuals 
# in HILO adaptation project
library(dplyr)
library(tidyr)
K = 2 # K = 4 number of genetic clusters/groups
# starting with pass1 analysis from 1st round of sequencing
admix<-read.table(paste0("../data/NGSadmix/merged_pass1_all_alloMaize4Low_16/K", K, "_pruned_all.qopt"))
colnames(admix) <- paste0("anc", 1:K) #c("anc1", "anc2")

# get metadata for individuals included in NGSadmix analysis
pass1_allo4Low <- read.table("../data/pass1_allo4Low_ids.txt", stringsAsFactors = F, 
                             header = T, sep = "\t")

# bams from pass1 matching NGSadmix output
#bams <- read.table("../data/pass1_allo4Low.all.list", stringsAsFactors = F) %>%
#  rename(bam = V1) %>%
#  mutate(n = as.numeric(gsub(".sort.dedup.baq.bam", "", 
#                             gsub("filtered_bam/hilo_", "", bam)))) %>%
#  left_join(., pass1, by = "n")
  
# define a threshold for high enough coverage to be included
#bamsPass <- bams[bams$est_coverage>=.2, ]
# join bams and admix by position (CAUTION - bam list order and admix results MUST MATCH!)
d <- bind_cols(pass1_allo4Low, admix)  %>%
  arrange(., popN) %>%
  arrange(., zea) %>%
  arrange(., symp_allo) %>%
  filter(., est_coverage > .25)

barplot(t(d[ , colnames(admix)]),col=rainbow(K),
        space=0,border=NA,
        xlab="Individuals",
        ylab="admixture")                

# now other subsets of data & SNP options
d1 <- read.table("../data/pass1_bam.all.list", stringsAsFactors = F) %>%
  rename(bam = V1) %>%
  mutate(n = as.numeric(gsub(".sort.dedup.baq.bam", "", 
                             gsub("filtered_bam/hilo_", "", bam)))) %>%
  bind_cols(., admix)  %>%
  arrange(., popN) %>%
  arrange(., zea) %>%
  arrange(., symp_allo)

admix_plot_prefix1 = "../plots/pass1_ngsadmix_SympVar_"
barplot(t(d1[ , c("V1", "V2")]),col=rainbow(K),
        space=0,border=NA,
        xlab="Individuals",
        ylab="admixture")
png(paste0(admix_plot_prefix1, "alloMex.png"), width = 800, height = 600)
barplot(t(d1[ d1$symp_allo == "allopatric" & d1$zea == "mexicana", c("V1", "V2")]),
        col=rainbow(K),
        space=0,border=NA,
        xlab="Individuals",
        ylab="admixture")
dev.off()
png(paste0(admix_plot_prefix1, "sympMaize.png"), width = 800, height = 600)
barplot(t(d1[ d1$symp_allo == "sympatric" & d1$zea == "maize", c("V1", "V2")]),
        col=rainbow(K),
        space=0,border=NA,
        xlab="Individuals",
        ylab="admixture")
dev.off()
png(paste0(admix_plot_prefix1, "sympMex.png"), width = 800, height = 600)
barplot(t(d1[ d1$symp_allo == "sympatric" & d1$zea == "mexicana", c("V1", "V2")]),
        col=rainbow(K),
        space=0,border=NA,
        xlab="Individuals",
        ylab="admixture")
dev.off()

# d1 subset is just allopatric mexicana and sympatric maize
infosub <- info[((info$zea == "mexicana" & info$symp_allo == "allopatric") | 
  (info$zea == "maize" & info$symp_allo == "sympatric")) & 
  info$est_coverage >= 0.05 & !is.na(info$est_coverage), ]
d1sub <- bind_cols(read.table(paste0("../data/NGSadmix/pass1/sympatricVar/SympMaizeAlloMex05/K", 
                                  K, "_pruned_all.qopt")), 
                infosub) %>%
                arrange(., popN) %>%
                arrange(., zea)
admix_plot_prefix = "../plots/pass1_ngsadmix_SympVar_SympMaizeAlloMex"
barplot(t(d1sub[ , c("V1", "V2")]),col=rainbow(K),
        space=0,border=NA,
        xlab="Individuals",
        ylab="admixture")
png(paste0(admix_plot_prefix, "alloMex.png"), width = 800, height = 600)
barplot(t(d1sub[ d1sub$symp_allo == "allopatric" & d1sub$zea == "mexicana", c("V1", "V2")]),
        col=rainbow(K),
        space=0,border=NA,
        xlab="Individuals",
        ylab="admixture")
dev.off()
png(paste0(admix_plot_prefix, "sympMaize.png"), width = 800, height = 600)
barplot(t(d1sub[ d1sub$symp_allo == "sympatric" & d1sub$zea == "maize", c("V1", "V2")]),
        col=rainbow(K),
        space=0,border=NA,
        xlab="Individuals",
        ylab="admixture")
dev.off()


# using original SNPs but just allopatric teosinte and sympatric maize
dsub <- bind_cols(read.table(paste0("../data/NGSadmix/pass1/allVar/SympMaizeAlloMex05/K", 
                                     K, "_pruned_all.qopt")), 
                   infosub) %>%
  arrange(., popN) %>%
  arrange(., zea)
admix_plot_prefix = "../plots/pass1_ngsadmix_allVar_SympMaizeAlloMex"
barplot(t(dsub[ , c("V1", "V2")]),col=rainbow(K),
        space=0,border=NA,
        xlab="Individuals",
        ylab="admixture")
png(paste0(admix_plot_prefix, "alloMex.png"), width = 800, height = 600)
barplot(t(dsub[ dsub$symp_allo == "allopatric" & dsub$zea == "mexicana", c("V1", "V2")]),
        col=rainbow(K),
        space=0,border=NA,
        xlab="Individuals",
        ylab="admixture")
dev.off()
png(paste0(admix_plot_prefix, "sympMaize.png"), width = 800, height = 600)
barplot(t(dsub[ dsub$symp_allo == "sympatric" & dsub$zea == "maize", c("V1", "V2")]),
        col=rainbow(K),
        space=0,border=NA,
        xlab="Individuals",
        ylab="admixture")
dev.off()




# allopatric populations
admix_plot_prefix = "../plots/pass1_ngsadmix_SNPs1_"
png(paste0(admix_plot_prefix, "alloMex.png"), width = 800, height = 600)
barplot(t(dPass[dPass$symp_allo == "allopatric", colnames(admix)]),col=rainbow(K),
        space=0,border=NA,
        xlab="Allopatric mexicana ind.'s",
        ylab="admixture")
dev.off()
# sympatric maize
png(paste0(admix_plot_prefix, "sympMaize.png"), width = 800, height = 600)
barplot(t(dPass[dPass$zea == "maize", colnames(admix)]),col=rainbow(K),
        space=0,border=NA,
        xlab="Sympatric maize ind.'s",
        ylab="admixture")
dev.off()
# sympatric mexicana
png(paste0(admix_plot_prefix, "sympMex.png"), width = 800, height = 600)
barplot(t(dPass[dPass$zea == "mexicana" & dPass$symp_allo == "sympatric", colnames(admix)]),col=rainbow(K),
        space=0,border=NA,
        xlab="Sympatric mexicana ind.'s",
        ylab="admixture")
dev.off()

# single pop visualize
barplot(t(dPass[dPass$popN == 22, colnames(admix)]),col=rainbow(K),
        space=0,border=NA,
        xlab="Individuals",
        ylab="admixture")
# ancestry 1 does not differentiate maize & mexicana well
tapply(dPass$anc1, dPass$popN, mean)
tapply(dPass$anc2, dPass$popN, mean)
# lets look at allele frequency estimates
freqs <- read.table("../data/NGSadmix/pass1/K2_pruned_all.fopt.gz") 
colnames(freqs) <- c("anc1", "anc2")
freqs$diff <- abs(freqs$anc1 - freqs$anc2)
freqs$mean <- (freqs$anc1 + freqs$anc2)/2
het1 <- mean(freqs$anc1 * (1-freqs$anc1))
het2 <- mean(freqs$anc2 * (1-freqs$anc2))
# a better total approximate would be to use the MAF of the entire sample estimated in ANGSD
# rather than the mean of two equal sized hypothetical pops
hetTot <- mean(freqs$mean * (1-freqs$mean))
# differentiation is fairly low between the two identified genetic clusters
Fst <- (hetTot - mean(het1, het2))/hetTot # [1] 0.03903411

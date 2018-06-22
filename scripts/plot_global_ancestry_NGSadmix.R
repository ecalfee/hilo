# calculation/plotting global ancestry for low-coverage maize/mexicana individuals 
# in HILO adaptation project
library(dplyr)
library(tidyr)
K = 2 # K = 4 number of genetic clusters/groups
# starting with pass1 analysis from 1st round of sequencing
admix<-read.table(paste0("../data/NGSadmix/pass1/K", K, "_pruned_all.qopt"))
colnames(admix) <- paste0("anc", 1:K) #c("anc1", "anc2")
# bams from pass1 matching NGSadmix output
bams <- read.table("../data/pass1_bam.all.list", stringsAsFactors = F)$V1
# get just HILO # from bam names
ids <- data.frame(n = as.numeric(gsub(".sort.dedup.baq.bam", "", 
                       gsub("filtered_bam/hilo_", "", bams))))

# overall sequencing/quality filtering 
qualUnsorted <- read.table("../data/filtered_bam/pass1.all.Nreads.sorted.csv", 
                   header = T, stringsAsFactors = F)
qualUnsorted$n <- as.numeric(gsub("hilo", "", qualUnsorted$ID))
qual <- arrange(qualUnsorted, n)
# histogram of coverage
png('../plots/pass1_hist_coverage_passQC.png', width = 800, height = 600)
hist(qual$est_coverage, breaks = 20, 
     main = "pass1 est. coverage post filtering", 
     xlab = "(# reads pass QC * 150bp) / 2.3Gb")
abline(v = median(qual$est_coverage), col = "blue")
abline(v = mean(qual$est_coverage), col = "red")
dev.off()

# number of individuals with higher coverage
sum(qual$est_coverage>.1)
sum(qual$est_coverage>.5)

# add in pass1 ind numbers (matches GL file, skips ind's with no seq data, not in 'qual'):
qual$pass1_ID <- paste0("Ind", 0:195)

meta <- read.table("../data/HILO_DAN_IDs_modEC.csv", stringsAsFactors = F)
colnames(meta) <- c("n", "ID", "fam")
meta <- separate(data = meta, col = fam, sep = "_", c("popN", "plate"), extra = "merge")
meta$popN <- as.numeric(meta$popN)
# label maize/mexicana and sympatric/allopatric using population number
meta$zea <- ifelse(meta$popN >= 100, "maize", "mexicana")
meta$symp_allo <- ifelse(meta$popN %in% c(20, 22, 33), "allopatric", "sympatric")

# print quality with IDs for all individuals in HILO_DAN_IDs_modEC.csv
info = left_join(meta,
                 qual[, c("num_read_pass", "est_coverage", "n", "pass1_ID")],
                 by = "n")
write.table(info, "../data/HILO_IDs_cov_pass1.csv", 
            quote = F, col.names = T, row.names = F)

# join meta data to list of included IDs
d <- left_join(ids, qual, by = "n") %>%
  left_join(., meta, by = "n") %>%
  bind_cols(admix, .) %>%
  arrange(., popN)
# define a threshold for high enough coverage to be included
dPass <- d[d$est_coverage>=.2,]
barplot(t(dPass[ , colnames(admix)]),col=rainbow(K),
        space=0,border=NA,
        xlab="Individuals",
        ylab="admixture")
# now other subsets of data & SNP options

d1 <- bind_cols(read.table(paste0("../data/NGSadmix/pass1/sympatricVar/allInd/K", 
                                K, "_pruned_all.qopt")), 
                info[!is.na(info$est_coverage), ]) %>%
                arrange(., popN) %>%
                arrange(., zea) %>%
                arrange(., symp_allo)
admix_plot_prefix = "../plots/pass1_ngsadmix_SympVar_"
barplot(t(d1[ , c("V1", "V2")]),col=rainbow(K),
        space=0,border=NA,
        xlab="Individuals",
        ylab="admixture")
png(paste0(admix_plot_prefix, "alloMex.png"), width = 800, height = 600)
barplot(t(d1[ d1$symp_allo == "allopatric" & d1$zea == "mexicana", c("V1", "V2")]),
        col=rainbow(K),
        space=0,border=NA,
        xlab="Individuals",
        ylab="admixture")
dev.off()
png(paste0(admix_plot_prefix, "sympMaize.png"), width = 800, height = 600)
barplot(t(d1[ d1$symp_allo == "sympatric" & d1$zea == "maize", c("V1", "V2")]),
        col=rainbow(K),
        space=0,border=NA,
        xlab="Individuals",
        ylab="admixture")
dev.off()
png(paste0(admix_plot_prefix, "sympMex.png"), width = 800, height = 600)
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

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
hist(qual$est_coverage)
# number of individuals with higher coverage
sum(qual$est_coverage>.1)
sum(qual$est_coverage>.5)

meta <- read.table("../data/HILO_DAN_IDs_modEC.csv", stringsAsFactors = F)
colnames(meta) <- c("n", "ID", "fam")
meta <- separate(data = meta, col = fam, sep = "_", c("popN", "plate"), extra = "merge")
meta$popN <- as.numeric(meta$popN)

# print quality with IDs for all individuals in HILO_DAN_IDs_modEC.csv
write.table(left_join(meta,
                      qual[, c("num_read_pass", "est_coverage", "n")],
                      by = "n"), "../data/HILO_IDs_cov_pass1.csv", 
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
# allopatric populations
barplot(t(dPass[dPass$popN %in% c(20, 22, 33), colnames(admix)]),col=rainbow(K),
        space=0,border=NA,
        xlab="Allopatric mexicana ind.'s",
        ylab="admixture")
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

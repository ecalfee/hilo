# calculation/plotting global ancestry for low-coverage maize/mexicana individuals 
# in HILO adaptation project
library(dplyr)
library(tidyr)
K = 2 # K = 4 number of genetic clusters/groups
# starting with pass1 analysis from 1st round of sequencing
file_prefix = paste0("../data/NGSadmix/merged_pass1_all_alloMaize4Low_16/K", K, "_pruned_all")
admix<-read.table(paste0(file_prefix, ".qopt"))
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

# lets look at allele frequency estimates for each ancestry
#freqs <- read.table("../data/NGSadmix/pass1/K2_pruned_all.fopt.gz") 
freqs <- read.table(paste0(file_prefix, ".fopt.gz"))
colnames(freqs) <- paste0("anc", 1:K)
freqs$diff <- abs(freqs$anc1 - freqs$anc2)
freqs$mean <- (freqs$anc1 + freqs$anc2)/2
het1 <- mean(freqs$anc1 * (1-freqs$anc1))
het2 <- mean(freqs$anc2 * (1-freqs$anc2))
# a better total approximate would be to use the MAF of the entire sample estimated in ANGSD
# rather than the mean of two equal sized hypothetical pops
hetTot <- mean(freqs$mean * (1-freqs$mean))
# differentiation is fairly low between the two identified genetic clusters
Fst <- (hetTot - mean(het1, het2))/hetTot # [1] 0.03903411
Fst # with allopatric maize included and using this new set of SNPs called including allopatric maize
# and filtered every 1000th SNP, Fst is much more reasonable [1] 0.1242707
  

# plot 'STRUCTURE-like' ancestry plots
png("../plots/NGSadmix_pass1.png", # saves plot as pin in ../plots/
    height = 5, width = 8, units = "in", res = 150)
d %>%
  select(., colnames(admix)) %>%
  t(.) %>%
  barplot(height = .,
          col=rainbow(K),
          space=0,
          border=NA,
          main = "all HILO populations",
          xlab="Individuals",
          ylab="admixture") 
mtext(paste("est. between-ancestry Fst:", round(Fst,2)))
dev.off()

# now plot again with individuals with very low coverage <0.05x filtered out      
png("../plots/NGSadmix_pass1_over_0.05x_coverage.png",
    height = 5, width = 8, units = "in", res = 150)
d %>%
  filter(., est_coverage >= .05) %>%
  select(., colnames(admix)) %>%
  t(.) %>%
  barplot(height = .,
          col=rainbow(K),
          space=0,
          border=NA,
          main = "all HILO populations",
          xlab="Individuals with est. coverage >= 0.05x",
          ylab="admixture") 
mtext(paste("est. between-ancestry Fst:", round(Fst,2)))
dev.off()

for (g in unique(d$group)){
  png(paste0("../plots/NGSadmix_pass1_", g, ".png"), # saves plot as pin in ../plots/
      height = 5, width = 8, units = "in", res = 150)
  d %>%
    filter(., group == g) %>%
    filter(., est_coverage >= .05) %>%
    select(., colnames(admix)) %>%
    t(.) %>%
    barplot(height = .,
            col=rainbow(K),
            space=0,
            border=NA,
            main = g,
            xlab="Individuals with est. coverage >= 0.05x",
            ylab="admixture")  
  dev.off()
}


# by location
p = d %>%
  filter(., est_coverage >= .05) %>%
  mutate(group = as.factor(paste(symp_allo, zea, sep = "_"))) %>%
  ggplot(., aes(x = group, y = anc1, color = group)) +
  geom_violin() + 
  geom_point(size=1, position = position_jitter(w=0.05)) +
  scale_colour_manual(values = c("yellow", "darkblue", "orange", "blue")) +
  theme_minimal() + 
  labs(y = "proportion 'mexicana-like' ancestry")
p + facet_wrap(~LOCALITY) # show plot
# save plot
ggsave("../plots/NGSadmix_pass1_proportion_mexicana-like_by_location.png", plot = p + facet_wrap(~LOCALITY), device = png(), 
       width = 12, height = 8, units = "in",
       dpi = 200)
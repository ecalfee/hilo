# this script evaluates pass1 population structure using various methods
# PCA of multiple SNP sets, Fst using all information from bam files (for a set of regions)
library(scales)
library(dplyr)
library(tidyr)
library(ggplot2)

# ID and population labels:
hilo <- read.table("../data/HILO_IDs_cov_pass1.csv", stringsAsFactors = F, header = T)
pass1 <- hilo[!is.na(hilo$num_read_pass),] # only include individuals with some pass1 data

# PCA
plot_PCA = function(file, name, jitter = F, ...){
  m <- as.matrix(read.table(file, stringsAsFactors = F, header = F))
  e <- eigen(m)
  colors = ifelse(pass1$zea=="maize", "orange", 
                  ifelse(pass1$symp_allo=="sympatric", "blue", "darkblue"))
  png(paste("../plots/pcangsd_bygroup", name, ".png"), # saves plot as pin in ../plots/
      height = 5, width = 8, units = "in", res = 150)
  if (jitter){ # plot with jitter
    plot(jitter(e$vectors[,1:2], .3), lwd=2, ylab="PC 2", xlab="PC 1",
         main=paste("PCA by group", name, "w/ jitter"),
         col = alpha(colors, 0.8),
         pch = 16, ...)
  } else{
    plot(e$vectors[,1:2], lwd=2, ylab="PC 2", xlab="PC 1",
         main=paste("PCA", name),
         col = alpha(colors, 0.8),
         pch = 16, ...)
  }
  legend("bottomleft", col = c("orange", "blue", "darkblue"), 
         legend = c("maize symp.", "mex. symp. ", "mex. allo."), 
         pch = 16, cex = .7)
  dev.off()
}
plot_PCA_byPop = function(file, name, jitter = F, ...){
  m <- as.matrix(read.table(file, stringsAsFactors = F, header = F))
  e <- eigen(m)
  largerainbow = rainbow(length(unique(pass1$popN)))
  colors = largerainbow[factor(pass1$popN)]
  points = ifelse(pass1$zea=="maize", 3, 
                  ifelse(pass1$symp_allo=="sympatric", 16, 1))
  png(paste("../plots/pcangsd_bypop", name, ".png"), # saves plot as pin in ../plots/
      height = 5, width = 8, units = "in", res = 150)
  plot(e$vectors[,1:2], lwd=2, ylab="PC 2", xlab="PC 1",
       main=paste("PCA", name),
       col = alpha(colors, 0.8),
       pch = points, ...)
  legend("bottomleft", 
         legend = c("maize symp.", "mex. symp. ", "mex. allo."), 
         pch = c(3, 16, 1), cex = .7)
  dev.off()
}

# PCAs -- keep in mind that of the 100K SNPs, several thousand were dropped in any pairwise Fst analysis
# (up to half in the smaller pop sizes) and many individuals have no data for PCA on some SNPs
plot_PCA(file = "../data/geno_lik/pass1/allVar/pruned_by1000.cov", name = "100K_SNPs")
plot_PCA(file = "../data/geno_lik/pass1/allVar/pruned_by1000.cov", name = "100K_SNPs_byCoverage", 
         cex = pass1$est_coverage)
plot_PCA(file = "../data/geno_lik/pass1/sympatricVar/pruned_all.cov", name = "SNPs_ascertained_sympatric_pops")

plot_PCA(file = "../data/geno_lik/pass1/pruned_all.cov", name = "SNPs_pruned_all")
plot_PCA(file = "../data/geno_lik/pass1/pruned_all.cov", name = "SNPs_pruned_all_jittered",
         jitter = T)
# 2 outliers here have very low coverage (disappear when I make point size ~ coverage)
plot_PCA(file = "../data/geno_lik/pass1/pruned_all.cov", name = "SNPs_pruned_all_jittered_byCoverage",
         jitter = T, cex = pass1$est_coverage)
# plotting first PCA again, by population, with no clear patterns
plot_PCA_byPop(file = "../data/geno_lik/pass1/allVar/pruned_by1000.cov", name = "100K_SNPs")

# plot pairwise Fst results here
fstG <- read.table("../data/FST/pass1/N1000.L100.regions/pairwise.group.fst.stats",
                  stringsAsFactors = F, header = T)
groupN <- c(1000, 2000, 3000)
groupName <- c("Maize_symp", "Mex_symp", "Mex_allo")
ggplot(data = fstG, 
       aes(x=pop1, y=pop2, fill=Fst_Hudson)) + 
  geom_tile()
mG <- matrix(0, 3, 3)
mG[1,2] <- fstG[fstG$pop1==1000 & fstG$pop2==2000, "Fst_Hudson"]
mG[1,3] <- fstG[fstG$pop1==1000 & fstG$pop2==3000, "Fst_Hudson"]
mG[2,3] <- fstG[fstG$pop1==2000 & fstG$pop2==3000, "Fst_Hudson"]
colnames(mG) <- groupName
rownames(mG) <- groupName
# Fst grid for maize & mexicana groups symp/allo
mG

# Below I get qualitatively similar conclusions from Fst using 1000 regions of 100 bp each 
# without needing to call variants
fst <- read.table("../data/FST/pass1/N1000.L100.regions/pairwise.fst.stats",
                  stringsAsFactors = F, header = T)
pops = unique(c(fst$pop1, fst$pop2))
fst$f1 <- as.numeric(factor(fst$pop1, levels = pops))
fst$f2 <- as.numeric(factor(fst$pop2, levels = pops))
m <- matrix(0, length(pops), length(pops))
for (i in 1:nrow(fst)){
  m[fst[i, "f1"], fst[i, "f2"]] <- fst[i, "Fst_Hudson"]
}
#corrplot(m, method = "shade")
ggplot(data = fst, 
       aes(x=f1, y=f2, fill=Fst_Hudson)) + 
  geom_tile() +
  ggtitle("all pops Fst")
ggplot(data = fst[fst$pop1 >= 100 & fst$pop2 >= 100, ], 
       aes(x=f1, y=f2, fill=Fst_Hudson)) + 
  geom_tile() +
  ggtitle("within maize Fst")
ggplot(data = fst[fst$pop1 %in% c(20,22,33) & fst$pop2 %in% c(20,22,33), ], 
       aes(x=f1, y=f2, fill=Fst_Hudson)) + 
  geom_tile() +
  ggtitle("within allo teosinte Fst")
ggplot(data = fst[fst$pop1 < 100 & !(fst$pop1 %in% c(20,22,33)) & fst$pop2 < 100 & !(fst$pop2 %in% c(20,22,33)), ],       
       aes(x=f1, y=f2, fill=Fst_Hudson)) + 
  geom_tile() +
  ggtitle("within symp. teosinte Fst")

# plot mapping statistics too; % mapped; % mapped > mapQ30 by mex/maize:
raw_metrics_no_header <- read.table("../data/filtered_bam/pass1.all.metrics.raw.Nreads", 
                          stringsAsFactors = F, header = F, skip = 1)
con = file("../data/filtered_bam/pass1.all.metrics.raw.Nreads", "r")
metrics_head <- strsplit(readLines(con, 1), split = "\t")[[1]]
close(con)
raw_metrics0 <- raw_metrics_no_header %>% 
  unite(., "LIB", V2, V3, sep= "_") # join 2 library_unknown columns
colnames(raw_metrics0) <- metrics_head # add header
raw_metrics <- raw_metrics0 %>%
  mutate(n_reads = 2*READ_PAIRS_EXAMINED + UNPAIRED_READS_EXAMINED) %>%
  mutate(prop_unmapped = UNMAPPED_READS/n_reads) %>%
  mutate(n_dedup = 2*(READ_PAIRS_EXAMINED - READ_PAIR_DUPLICATES) +
           UNPAIRED_READS_EXAMINED - UNPAIRED_READ_DUPLICATES) %>%
  mutate(n = as.integer(substr(StudyID, 5, 10))) %>%
  left_join(., hilo, by = "n") %>%
  mutate(n_filtOutQ30 = n_dedup-num_read_pass) %>%
  mutate(prop_filtOutQ30 = n_filtOutQ30/n_dedup) %>%
  mutate(prop_filtTot = (n_reads - num_read_pass)/n_reads) %>%
  mutate(n_dropped = n_reads-n_dedup+n_filtOutQ30)


# plots for coverage and filtering effects
png(paste("../plots/pass1_total_reads_dropped_by_group.png"), # saves plot as pin in ../plots/
    height = 6, width = 6, units = "in", res = 150)
qplot(prop_filtTot, data=raw_metrics, geom="density", fill=paste(zea, symp_allo, sep = "_"), alpha=I(.5), 
      main="Prop. reads filtered total by group", xlab="Prop. reads dropped", 
      ylab="Density") +
  theme(legend.title=element_blank())
dev.off()
# same but in scatterplot across sequencing coverage
ggplot(data = raw_metrics) + geom_point(aes(x = n_reads, 
                                            y = prop_filtTot,
                                            color = paste(zea, symp_allo, sep = "_"))) +
  theme(legend.title=element_blank())
png(paste("../plots/pass1_reads_unmapped_by_group.png"), # saves plot as pin in ../plots/
    height = 6, width = 6, units = "in", res = 150)
qplot(prop_unmapped, data=raw_metrics, geom="density", fill=paste(zea, symp_allo, sep = "_"), alpha=I(.5), 
      main="Prop. unmapped reads by group", xlab="Prop. reads unmapped", 
      ylab="Density") +
  theme(legend.title=element_blank())
dev.off()
ggplot(data = raw_metrics) + geom_point(aes(x = n_reads, 
                                            y = prop_unmapped,
                                            color = paste(zea, symp_allo, sep = "_"))) +
  theme(legend.title=element_blank())
# of reads that mapped and are not duplicates -- what proportion are mapQ <30?
png(paste("../plots/pass1_prop_reads_mapped_but_belowQ30_by_group.png"), # saves plot as pin in ../plots/
    height = 6, width = 6, units = "in", res = 150)
qplot(prop_filtOutQ30, data=raw_metrics, geom="density", fill=paste(zea, symp_allo, sep = "_"), alpha=I(.5), 
      main="Prop. reads filtered by mapQ<30 by group", xlab="Prop. reads < Q30", 
      ylab="Density") +
  theme(legend.title=element_blank())
dev.off()
# strangely, teosinte allopatric has higher percent duplications
png(paste("../plots/pass1_prop_duplication_by_group.png"), # saves plot as pin in ../plots/
    height = 6, width = 6, units = "in", res = 150)
qplot(PERCENT_DUPLICATION, data=raw_metrics, geom="density", fill=paste(zea, symp_allo, sep = "_"), alpha=I(.5), 
      main="Prop. reads duplications by group", xlab="Prop. reads duplicates", 
      ylab="Density") +
  theme(legend.title=element_blank())
dev.off()

# unclear what the effect of BAQ (adjusted base quality) filtering is .. could try to count for a set of SNPs
# Jeff thinks a larger issue is getting total coverage across individuals, groups and entire sample for SNPs included (e.g. in PCA)

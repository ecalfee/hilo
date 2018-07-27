# this script evaluates pass1 population structure using various methods
# PCA of multiple SNP sets, Fst using all information from bam files (for a set of regions)
library(scales)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

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
fstG <- read.table("../data/SFS/pass1/N1000.L100.regions/pairwise.group.fst.stats",
                  stringsAsFactors = F, header = T)
groupN <- c(1000, 2000, 3000, 5000)
groupName <- c("Maize_symp", "Mex_symp", "Mex_allo", "Maize_allo_lowland")
#ggplot(data = fstG, 
#       aes(x=pop1, y=pop2, fill=Fst_Hudson)) + 
#  geom_tile()
mG <- matrix(0, length(groupN), length(groupN))
for (i in fstG[, "pop1"]){
  for (j in fstG[fstG$pop1 == i, "pop2"]){
    mG[which(groupN == i), which(groupN == j)] <- fstG[fstG$pop1 == i & fstG$pop2 == j, "Fst_Hudson"]
    mG[which(groupN == j), which(groupN == i)] <- fstG[fstG$pop1 == i & fstG$pop2 == j, "Fst_Hudson"]
  }
}
colnames(mG) <- groupName
rownames(mG) <- groupName
# Fst grid for maize & mexicana groups symp/allo
mG <- mG[c(4,1,2,3), c(4,1,2,3)] # reorder to go allopatric maize to allopatric mexicana

# plot Fst nicely
data_mG <- reshape2::melt(mG[, ]) %>%
  rename(., Group1=Var1, Group2=Var2)

p <- ggplot(data =  data_mG, aes(x = Group1, y = Group2)) +
  geom_tile(aes(fill = value), colour = "white") +
  geom_text(aes(label = sprintf("%1.3f",value)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ggtitle("Fst (Hudson's) from 1000 regions, each 100bp") + 
  theme(axis.title.x = element_blank(), 
       axis.title.y = element_blank())
p
ggsave("../plots/Fst_maize_mex_100000bp.png", plot = p, device = png(), 
       width = 10, height = 8, units = "in",
       dpi = 100)


# Below I get qualitatively similar conclusions from Fst using 1000 regions of 100 bp each 
# without needing to call variants. All subpopulations plotted below.
fst <- read.table("../data/SFS/pass1/N1000.L100.regions/pairwise.fst.stats",
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


# this script evaluates pass1 population structure using various methods
# PCA of multiple SNP sets, Fst using all information from bam files (for a set of regions)
library(scales)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

# ID and population labels:
hilo <- read.table("../data/HILO_IDs_cov_pass1.csv", stringsAsFactors = F, header = T)
# add in metadata
# germplasm data from JRI for all projects ("riplasm"). 
# RIMMA is hilo maize and RIMME is hilo mex. This just gives me metadata for each population..(e.g. lat/long)
riplasm <- read.csv("../data/riplasm/riplasm.csv", stringsAsFactors = F, header = T) %>%
  mutate(., prefix = substr(RI_ACCESSION, 1, 5)) %>%
  mutate(., ind = as.integer(substr(RI_ACCESSION, 6, 100))) %>%
  dplyr::filter(., prefix %in% c("RIMMA", "RIMME"))
hilo_meta <- cbind(hilo, do.call(rbind,
                                 apply(hilo, 1, function(row) riplasm[riplasm$RI_ACCESSION ==  paste0(
                                   ifelse(row["zea"] == "mexicana", "RIMME", "RIMMA"), 
                                   sprintf("%04d", as.integer(row["popN"]))), 
                                   c("RI_ACCESSION", "GEOCTY", "LOCALITY")])))

# metadata for 16 individuals from 4 lowland populations
allo4Low <- data.frame(popN=c(0,0,0,0,-10,-10,-10,-10,-20,-20,-20,-20,-30,-30,-30,-30),
                     zea = rep("maize", 16),
                     symp_allo = rep("allopatric", 16),
                     ID = paste0("4Low", c(1:4,11:14,21:24,31:34)),
                     est_coverage = rep(2, 16), # underestimate of coverage (but ok for now)
                     RI_ACCESSION = NA,
                     GEOCTY = "Mexico",
                     LOCALITY = "Lowland_4pops",
                     stringsAsFactors = F)
pass1 <- hilo_meta[!is.na(hilo_meta$num_read_pass),] # only include individuals with some pass1 data
pass1_allo4Low <- select(pass1, c("ID", "popN","zea", "symp_allo", "est_coverage", "RI_ACCESSION", "GEOCTY", "LOCALITY")) %>%
  bind_rows(., allo4Low)

metrics <- read.table("../data/filtered_bam/pass1.all.metrics.calcs", stringsAsFactors = F, header = T)

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
# more restrictive SNP pruning that takes only variants wiht >100 ind's with data (more than half)
plot_PCA(file = "../data/geno_lik/pass1/pruned_all.cov", name = "SNPs_pruned_all")
plot_PCA(file = "../data/geno_lik/pass1/pruned_all.cov", name = "SNPs_pruned_all_jittered",
         jitter = T)
# 2 outliers here have very low coverage (disappear when I make point size ~ coverage)
plot_PCA(file = "../data/geno_lik/pass1/pruned_all.cov", name = "SNPs_pruned_all_jittered_byCoverage",
         jitter = T, cex = pass1$est_coverage)
# plotting first PCA again, by population, with no clear patterns
plot_PCA_byPop(file = "../data/geno_lik/pass1/allVar/pruned_by1000.cov", name = "100K_SNPs")

# more PCA's
# just first 40 individuals, then first 98 (plates 1 & 2) look more as expected but don't have allopatric teosinte ind's:
start_n = 1
end_n = 98
#end_n = 196
m <- as.matrix(read.table(file =  "../data/geno_lik/pass1/allVar/pruned_by1000.cov", 
                          stringsAsFactors = F, header = F))[start_n:end_n, start_n:end_n]
e <- eigen(m)
colors = ifelse(pass1$zea[start_n:end_n]=="maize", "orange", 
                ifelse(pass1$symp_allo[start_n:end_n]=="sympatric", "blue", "darkblue"))
plot(e$vectors[,1:2], lwd=2, ylab="PC 2", xlab="PC 1",
     main=paste("PCA", "of HILO id's", start_n, "to", end_n),
     col = alpha(colors, 0.8),
     pch = 16,
     cex = pass1$est_coverage[start_n:end_n])
# additional PC's
plot(e$vectors[,3:4], lwd=2, ylab="PC 4", xlab="PC 3",
     main=paste("PCA", "of HILO id's", start_n, "to", end_n),
     col = alpha(colors, 0.8),
     pch = 16,
     cex = pass1$est_coverage[start_n:end_n])
plot(e$vectors[,6:7], lwd=2, ylab="PC 7", xlab="PC 6",
     main=paste("PCA", "of HILO id's", start_n, "to", end_n),
     col = alpha(colors, 0.8),
     pch = 16,
     cex = pass1$est_coverage[start_n:end_n])
  legend("bottomleft", col = c("orange", "blue", "darkblue"), 
         legend = c("maize symp.", "mex. symp. ", "mex. allo."), 
         pch = 16, cex = .7)

# visualize pruned SNPs with over 100 ind's with data without jitter (Removing outliers)
  # this is more restrictive SNP pruning
  min_cov = .2
  m1 <- as.matrix(read.table(file =  "../data/geno_lik/pass1/pruned_all.cov", 
                            stringsAsFactors = F, header = F))[pass1$est_coverage >= min_cov, pass1$est_coverage >= min_cov]
  e1 <- eigen(m1)
  colors1 = ifelse(pass1$zea[pass1$est_coverage >= min_cov]=="maize", "orange", 
                  ifelse(pass1$symp_allo[pass1$est_coverage >= min_cov]=="sympatric", "blue", "darkblue"))
  # If I include low coverage individuals, all main PC's are just tagging outliers
  # otherwise, it gets a similar spread to the other PCA (e.g. if including only > .2 ind coverage)
  plot(e1$vectors[,1:2], lwd=2, ylab="PC 2", xlab="PC 1",
       main=paste("PCA", "of HILO id's w/ min. cov.=", min_cov),
       col = alpha(colors1, 0.8),
       pch = 16,
       cex = 5*pass1$est_coverage[pass1$est_coverage >= min_cov])
  plot(e1$vectors[,3:4], lwd=2, ylab="PC 4", xlab="PC 3",
       main=paste("PCA", "of HILO id's w/ min. cov.=", min_cov),
       col = alpha(colors1, 0.8),
       pch = 16,
       cex = 5*pass1$est_coverage[pass1$est_coverage >= min_cov])
  plot(e1$vectors[,5:6], lwd=2, ylab="PC 6", xlab="PC 5",
       main=paste("PCA", "of HILO id's w/ min. cov.=", min_cov),
       col = alpha(colors1, 0.8),
       pch = 16,
       cex = 5*pass1$est_coverage[pass1$est_coverage >= min_cov])
  plot(e1$vectors[,11:12], lwd=2, ylab="PC 12", xlab="PC 11",
       main=paste("PCA", "of HILO id's w/ min. cov.=", min_cov),
       col = alpha(colors1, 0.8),
       pch = 16,
       cex = 5*pass1$est_coverage[pass1$est_coverage >= min_cov])
  legend("bottomleft", col = c("orange", "blue", "darkblue"), 
         legend = c("maize symp.", "mex. symp. ", "mex. allo."), 
         pch = 16, cex = .7)  
  

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


# ADDING IN ALLOPATRIC MAIZE - 16 individuals from 4 lowland populations
#m2 <- as.matrix(read.table("../data/geno_lik/merged_pass1_all_alloMaize4Low_16/allVar/whole_genome_pruned_every_1000.partial.cov",
#                          stringsAsFactors = F, header = F))
m2 <- as.matrix(read.table("../data/geno_lik/merged_pass1_all_alloMaize4Low_16/allVar/whole_genome_pruned_every_1000.cov",
                           stringsAsFactors = F, header = F))
e2 <- eigen(m2)
colors2 = ifelse(pass1_allo4Low$zea=="maize", ifelse(pass1_allo4Low$symp_allo=="sympatric", "orange", "yellow"), 
                ifelse(pass1_allo4Low$symp_allo=="sympatric", "blue", "darkblue"))
png(paste("../plots/pcangsd_bygroup", "with_16_lowland_maize", ".png"), # saves plot as pin in ../plots/
    height = 5, width = 8, units = "in", res = 300)
plot(e2$vectors[,1:2], lwd=2, ylab=paste("PC 2 -",e2$values[2], "%"), xlab=paste("PC 1 -", e2$values[1], "%"),
       main=paste("PCA adding lowland maize"),
       col = alpha(colors2, 0.8),
       pch = 16)
legend("topright", col = c("yellow", "orange", "blue", "darkblue"), 
       legend = c("maize allo.", "maize symp.", "mex. symp. ", "mex. allo."), 
       pch = 16, cex = .7)
dev.off()

# point size by sequencing amount
png(paste("../plots/pcangsd_bygroup_seqDepth", "with_16_lowland_maize", ".png"), # saves plot as pin in ../plots/
    height = 5, width = 8, units = "in", res = 300)
plot(e2$vectors[,1:2], lwd=2, ylab=paste("PC 2 -",e2$values[2], "%"), xlab=paste("PC 1 -", e2$values[1], "%"),
     main=paste("PCA adding lowland maize; size = depth"),
     col = alpha(colors2, 0.8),
     pch = 16,
     #cex = c(mean_depthF, rep(2, 16))*2)
     cex = pass1_allo4Low$est_coverage*2)
legend("topright", col = c("yellow", "orange", "blue", "darkblue"), 
       legend = c("maize allo.", "maize symp.", "mex. symp. ", "mex. allo."), 
       pch = 16, 
       cex = 0.7)
dev.off()

# plot with allopatric maize by group
png(paste("../plots/pcangsd_bypop", "with_16_lowland_maize", ".png"), # saves plot as pin in ../plots/
    height = 5, width = 8, units = "in", res = 300)
largerainbow2 = rainbow(length(unique(pass1_allo4Low$popN)))
colorsrainbow2 = largerainbow2[factor(pass1_allo4Low$popN)]
points2 = ifelse(pass1_allo4Low$zea=="maize", ifelse(pass1_allo4Low$symp_allo=="sympatric", 3, 9), 
                 ifelse(pass1_allo4Low$symp_allo=="sympatric", 16, 1))
plot(e2$vectors[,1:2], lwd=2, ylab=paste("PC 2 -",e2$values[2], "%"), xlab=paste("PC 1 -", e2$values[1], "%"),
     main=paste("PCA adding lowland maize"),
     col = alpha(colorsrainbow2, 0.8),
     pch = points2)
legend("topright",
       legend = c("maize allo.", "maize symp.", "mex. symp. ", "mex. allo."), 
       pch = c(9, 3, 16, 1), cex = .7)
dev.off()

# plot just first 40 sequenced individuals:
sub = c(1:40,197:212)
e2sub <- eigen(m2[sub, sub])
png(paste("../plots/pcangsd_bygroup_hilo1-40_seqDepth", "with_16_lowland_maize", ".png"), # saves plot as pin in ../plots/
    height = 5, width = 8, units = "in", res = 150)
plot(e2sub$vectors[,1:2], lwd=2, ylab=paste("PC 2 -",e2sub$values[2], "%"), xlab=paste("PC 1 -", e2sub$values[1], "%"),
     main=paste("PCA Hilo 1-40 only adding lowland maize; size = depth"),
     col = alpha(colors2[sub], 0.8),
     pch = 16,
     cex = pass1_allo4Low$est_coverage[sub]*2)
legend("topright", col = c("yellow", "orange", "blue", "darkblue"), 
       legend = c("maize allo.", "maize symp.", "mex. symp. ", "mex. allo."), 
       pch = 16, 
       cex = 0.7)
dev.off()

# write new file with metadata and PC's 1 and 2 for hilo individuals and 16 lowland allopatric maize too
pc2 = e2$vectors[,1:2]
colnames(pc2) = c("PC1", "PC2")
pass1_allo4Low_wPC <- cbind(pass1_allo4Low, pc2)
write.table(pass1_allo4Low_wPC, "../data/geno_lik/merged_pass1_all_alloMaize4Low_16/allVar/metadata_and_whole_genome_pruned_every_1000_2PCs.txt", sep = "\t",
            col.names = T, row.names = F, quote = F)

pass1_allo4Low_wPC_wMetrics = left_join(pass1_allo4Low_wPC, metrics, 
             by = c("ID", "popN", "zea", "symp_allo", "est_coverage"))

pass1_allo4Low_wPC_wMetrics %>%
  ggplot(., aes(PC1, prop_unmapped)) +
  geom_point(aes(color = paste(symp_allo, zea, sep = "_"),
                 size = depthRegions)) +
                 #size = est_coverage)) +
  scale_colour_manual(values = c("yellow", "darkblue", "orange", "blue")) +
  labs(color = "Group", size = "Est. coverage", y = "proportion reads not mapped")
ggsave("../plots/PC1_by_unmapped.png", device = png(), 
       width = 8, height = 6, units = "in",
       dpi = 200)

pass1_allo4Low_wPC_wMetrics %>%
  ggplot(., aes(PC1, prop_filtTot)) +
  geom_point(aes(color = paste(symp_allo, zea, sep = "_"), 
                 size = depthRegions)) +
  #size = est_coverage)) +
  scale_colour_manual(values = c("yellow", "darkblue", "orange", "blue")) +
  labs(color = "Group", size = "Est. coverage", y = "proportion reads filtered total")
ggsave("../plots/PC1_by_filtered_total.png", device = png(), 
       width = 8, height = 6, units = "in",
       dpi = 200)

pass1_allo4Low_wPC_wMetrics %>%
  ggplot(., aes(PC1, prop_filtOutQ30)) +
  geom_point(aes(color = paste(symp_allo, zea, sep = "_"), 
                 size = depthRegions)) +
  #size = est_coverage)) +
  scale_colour_manual(values = c("yellow", "darkblue", "orange", "blue")) +
  labs(color = "Group", size = "Est. coverage", y = "proportion reads mapped but mapQ<30")
ggsave("../plots/PC1_by_mapQ30.png", device = png(), 
       width = 8, height = 6, units = "in",
       dpi = 200)

# plot by location
p <- ggplot(pass1_allo4Low_wPC_wMetrics, aes(PC1, PC2)) + 
  geom_point(aes(color = paste(symp_allo, zea, sep = "_"), 
                 size = est_coverage)) +
  scale_colour_manual(values = c("yellow", "darkblue", "orange", "blue")) +
  labs(color = "Group", size = "Est. coverage")

# Use vars() to supply faceting variables:
p + facet_wrap(~LOCALITY)
ggsave("../plots/PCA_facet_by_population.png", plot = p + facet_wrap(~LOCALITY), device = png(), 
       width = 12, height = 8, units = "in",
       dpi = 200)







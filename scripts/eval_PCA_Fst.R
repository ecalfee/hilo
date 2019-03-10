# this script evaluates pass1 population structure using various methods
# PCA of multiple SNP sets, Fst using all information from bam files (for a set of regions)
library(scales)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

# ID and population labels from addMetadata2HiloIDs.R:
pass1 <- read.table("../data/pass1_ids.txt", stringsAsFactors = F, header = T, sep = "\t")

# with allopatric maize too
pass1_allo4Low <- read.table("../data/pass1_allo4Low_ids.txt", stringsAsFactors = F, header = T, sep = "\t")

# add metrics from plot_mapping_metrics.R
metrics <- read.table("../data/filtered_bam/pass1.all.metrics.calcs", stringsAsFactors = F, header = T,
                      sep = "\t")

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
#start_n = 161
start_n = 1
end_n = 196
#end_n = 196
#end_n = 160
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
     cex = 3*pass1$est_coverage[start_n:end_n])
# additional PC's
plot(e$vectors[,3:4], lwd=2, ylab="PC 4", xlab="PC 3",
     main=paste("PCA", "of HILO id's", start_n, "to", end_n),
     col = alpha(colors, 0.8),
     pch = 16,
     cex = 3*pass1$est_coverage[start_n:end_n])
plot(e$vectors[,6:7], lwd=2, ylab="PC 7", xlab="PC 6",
     main=paste("PCA", "of HILO id's", start_n, "to", end_n),
     col = alpha(colors, 0.8),
     pch = 16,
     cex = 3*pass1$est_coverage[start_n:end_n])
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
       cex = 3*pass1$est_coverage[pass1$est_coverage >= min_cov])
  plot(e1$vectors[,3:4], lwd=2, ylab="PC 4", xlab="PC 3",
       main=paste("PCA", "of HILO id's w/ min. cov.=", min_cov),
       col = alpha(colors1, 0.8),
       pch = 16,
       cex = 3*pass1$est_coverage[pass1$est_coverage >= min_cov])
  plot(e1$vectors[,5:6], lwd=2, ylab="PC 6", xlab="PC 5",
       main=paste("PCA", "of HILO id's w/ min. cov.=", min_cov),
       col = alpha(colors1, 0.8),
       pch = 16,
       cex = 3*pass1$est_coverage[pass1$est_coverage >= min_cov])
  plot(e1$vectors[,11:12], lwd=2, ylab="PC 12", xlab="PC 11",
       main=paste("PCA", "of HILO id's w/ min. cov.=", min_cov),
       col = alpha(colors1, 0.8),
       pch = 16,
       cex = 3*pass1$est_coverage[pass1$est_coverage >= min_cov])
  legend("bottomleft", col = c("orange", "blue", "darkblue"), 
         legend = c("maize symp.", "mex. symp. ", "mex. allo."), 
         pch = 16, cex = .7)  
  

# plot pairwise Fst results here
fstG <- read.table("../data/SAF/pass1/N1000.L100.regions/pairwise.group.fst.stats",
                  stringsAsFactors = F, header = F)
colnames(fstG) <- c("pop1", "pop2", "Fst_other", "Fst_Hudson")
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
ggsave("../plots/Fst_maize_mex_100000bp.png", plot = p, device = "png", 
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
# write data frame to file to share
write.table(pass1_allo4Low_wPC, "../data/geno_lik/merged_pass1_all_alloMaize4Low_16/allVar/metadata_and_whole_genome_pruned_every_1000_2PCs_corrected8.8.18.txt", sep = "\t",
            col.names = T, row.names = F, quote = F)

# get list of populations & familys (unique ID's individuals)
# grown in second greenhouse batch together
grow2 = read.csv("../data/pre_label_fix/NewAccessionsForSequencing_fromDan.csv", 
                 header = T, 
                 stringsAsFactors = F)[, c("Population", "Family")]%>%
  mutate(., grow = "greenhouse2") %>%
  rename(., popN = Population) %>%
  rename(., family = Family)

pass1_allo4Low_wPC_grow2 <- left_join(pass1_allo4Low_wPC, grow2,
                                      by = c("popN", "family")) %>%
  mutate(grow = ifelse(is.na(grow), "other", grow)) %>%
  mutate(., cluster = ifelse(PC1 < .025 & PC1 > -.01, "unclear",
                             ifelse((zea == "mexicana" & PC1 < 0) | 
                                      (zea == "maize" & PC1 > 0),
                                    "withOther", "withSame")))
# write data frame to file to share which individuals cluster differently
write.table(pass1_allo4Low_wPC_grow2, "../data/geno_lik/merged_pass1_all_alloMaize4Low_16/allVar/metadata_and_whole_genome_pruned_every_1000_2PCs_clusterGroup_corrected8.8.18.txt", sep = "\t",
            col.names = T, row.names = F, quote = F)


# plot highlighting individuals that are identified as clustering with 'same' or 'other' group
# or are intermediate 'unclear'
pass1_allo4Low_wPC_grow2 %>%
  ggplot(., aes(PC1, PC2)) +
  geom_point(aes(color = paste(symp_allo, zea, sep = "_"),
                 size = est_coverage, shape = cluster)) +
  scale_colour_manual(values = c("yellow", "darkblue", "orange", "blue")) +
  labs(color = "Group", size = "Est. coverage")

# plot with highlighting individuals from greenhouse2
pass1_allo4Low_wPC_grow2 %>%
  ggplot(., aes(PC1, PC2)) +
  geom_point(aes(color = paste(symp_allo, zea, sep = "_"),
                 size = est_coverage)) +
  scale_colour_manual(values = c("yellow", "darkblue", "orange", "blue")) +
  labs(color = "Group", size = "Est. coverage") +
  facet_wrap(~as.factor(grow))
ggsave("../plots/PCA_by_greenhouse_group.png", device = "png", 
       width = 12, height = 6, units = "in",
       dpi = 200)
# plotting again, this time so it's easier to see the one outlier
pass1_allo4Low_wPC_grow2 %>%
  ggplot(., aes(PC1, PC2)) +
  geom_point(aes(color = grow,
                 size = est_coverage,
                 alpha = .6)) +
  labs(color = "Grow Group", size = "Est. coverage") +
  facet_wrap(~zea)
ggsave("../plots/PCA_by_greenhouse_group2.png", device = "png", 
       width = 12, height = 6, units = "in",
       dpi = 200)

# after getting new data from anne, re-plotting PCA
pass1_allo4Low_wPC_grow2 %>%
  ggplot(., aes(PC1, PC2)) +
  geom_point(aes(color = paste(symp_allo, zea, sep = "_"),
                 size = est_coverage)) +
  scale_colour_manual(values = c("yellow", "darkblue", "orange", "blue")) +
  labs(color = "Group", size = "Est. coverage")
ggsave("../plots/PCA_newIDs_all.png", device = "png", 
       width = 12, height = 6, units = "in",
       dpi = 200)
# exclude individuals with very low coverage < 0.5x estimated from # reads passing sequencing fiilters
pass1_allo4Low_wPC_grow2 %>%
  filter(est_coverage > .05) %>%
  ggplot(., aes(PC1, PC2)) +
  geom_point(aes(color = paste(symp_allo, zea, sep = "_"),
                 size = est_coverage)) +
  scale_colour_manual(values = c("yellow", "darkblue", "orange", "blue")) +
  labs(color = "Group", size = "Est. coverage")
ggsave("../plots/PCA_newIDs_all_over_0.05x_coverage.png", device = "png", 
       width = 12, height = 6, units = "in",
       dpi = 200)


# separating out individuals 143-160 with new ID's from Anne
pass1_allo4Low_wPC_grow2 %>%
  mutate(id_batch = ifelse(grow == "greenhouse2" & n <= 160, "hilo143-160", "hilo <142 or >160")) %>%
  ggplot(., aes(PC1, PC2)) +
  geom_point(aes(color = paste(symp_allo, zea, sep = "_"),
                 size = est_coverage)) +
  scale_colour_manual(values = c("yellow", "darkblue", "orange", "blue")) +
  labs(color = "Group", size = "Est. coverage") +
  facet_wrap(~as.factor(id_batch))
ggsave("../plots/PCA_newIDs_sep143-160.png", device = "png", 
       width = 12, height = 6, units = "in",
       dpi = 200)
# same plot but with individuals given different shapes for different populations
for (i in c("first142", "hilo <143 or >160", "hilo > 160")){
  pass1_allo4Low_wPC_grow2 %>%
    mutate(id_batch = ifelse(n < 143, "first142", ifelse(grow == "greenhouse2" & n <= 160, 
                                                         "hilo143-160", "hilo <143 or >160"))) %>%
    filter(., id_batch == i) %>%
    ggplot(., aes(PC1, PC2)) +
    geom_point(aes(shape = as.factor(id_batch),
                   size = est_coverage,
                   color = as.factor(popN))) +
    labs(color = "Population N", size = "Est. coverage", shape = "ID set") +
    facet_wrap(~as.factor(paste(symp_allo, zea, sep = "_")))
}

ggsave("../plots/PCA_newIDs_sep143-160_by population.png", device = "png", 
       width = 12, height = 6, units = "in",
       dpi = 200)


pass1_allo4Low_wPC_wMetrics = left_join(pass1_allo4Low_wPC, metrics, 
             by = c("ID", "popN", "zea", "symp_allo", "est_coverage"))

pass1_allo4Low_wPC_wMetrics %>%
  ggplot(., aes(PC1, prop_unmapped)) +
  geom_point(aes(color = paste(symp_allo, zea, sep = "_"),
                 size = depthRegions)) +
                 #size = est_coverage)) +
  scale_colour_manual(values = c("yellow", "darkblue", "orange", "blue")) +
  labs(color = "Group", size = "Est. coverage", y = "proportion reads not mapped")
ggsave("../plots/PC1_by_unmapped.png", device = "png", 
       width = 8, height = 6, units = "in",
       dpi = 200)

pass1_allo4Low_wPC_wMetrics %>%
  ggplot(., aes(PC1, prop_filtTot)) +
  geom_point(aes(color = paste(symp_allo, zea, sep = "_"), 
                 size = depthRegions)) +
  #size = est_coverage)) +
  scale_colour_manual(values = c("yellow", "darkblue", "orange", "blue")) +
  labs(color = "Group", size = "Est. coverage", y = "proportion reads filtered total")
ggsave("../plots/PC1_by_filtered_total.png", device = "png", 
       width = 8, height = 6, units = "in",
       dpi = 200)

pass1_allo4Low_wPC_wMetrics %>%
  ggplot(., aes(PC1, prop_filtOutQ30)) +
  geom_point(aes(color = paste(symp_allo, zea, sep = "_"), 
                 size = depthRegions)) +
  #size = est_coverage)) +
  scale_colour_manual(values = c("yellow", "darkblue", "orange", "blue")) +
  labs(color = "Group", size = "Est. coverage", y = "proportion reads mapped but mapQ<30")
ggsave("../plots/PC1_by_mapQ30.png", device = "png", 
       width = 8, height = 6, units = "in",
       dpi = 200)

# plot by location
p <- pass1_allo4Low_wPC_grow2 %>%
  filter(n <143 | n > 160) %>%
  ggplot(., aes(PC1, PC2)) + 
  geom_point(aes(color = paste(symp_allo, zea, sep = "_"), 
                 size = est_coverage)) +
  scale_colour_manual(values = c("darkblue", "orange", "blue")) +
  labs(color = "Group", size = "Est. coverage")

# Use vars() to supply faceting variables:
p + facet_wrap(~LOCALITY)
ggsave("../plots/PCA_facet_by_population.png", plot = p + facet_wrap(~LOCALITY), device = "png", 
       width = 12, height = 8, units = "in",
       dpi = 200)

# plot PCA for outliers
m <- as.matrix(read.table("../data/outliers/chr4/inv4m/PCA.cov", stringsAsFactors = F, header = F))
e <- eigen(m)
colors = ifelse(pass1$zea=="maize", "orange", 
                ifelse(pass1$symp_allo=="sympatric", "blue", "darkblue"))

# load inv4m_anc from plotLocalAncTracts.R for labels
pca = e$vectors[, 1:2]
colnames(pca) = paste0("PC", 1:2)
pca2 = cbind(pass1_allo4Low, pca) %>%
  left_join(., inv4m_anc, by = c("ID", "n"="hilo_id", "zea"))
# there are some NAs because they are from pops without ancestry calls
pca2$hap_group[is.na(pca2$hap_group)] <- paste(pca2$zea, pca2$symp_allo, sep = "_")[is.na(pca2$hap_group)]
pca2 %>%
  ggplot(., aes(x = PC1, y = PC2)) +
  geom_point(aes(color = hap_group))
pca2 %>%
  #filter(!hap_group == "maize_sympatric" & !hap_group == "mexicana_sympatric") %>%
  #filter(!LOCALITY == "Puerta Encantada") %>%
  ggplot(., aes(x = PC1, y = PC2)) +
  geom_point(aes(color = anc, 
                 #size = top_anc_n,
                 shape = paste(zea,symp_allo)),
             alpha = .8) + 
  labs(title = "PCA of haplotypes at inv4m",
       x = paste("PC 1 -", round(e$values[1]), "%"),
       y = paste("PC 2 -", round(e$values[2]), "%"))
ggsave("../plots/PCA_inv4m.png", device = "png", 
       width = 12, height = 10, units = "in",
       dpi = 200)
# do the mexicana haplotypes appear to cluster at all by location?
# no, not really. I see this when I just plot mex.mex haplotypes
pca2 %>%
  filter(anc == "mex.mex" | hap_group == "mexicana_allopatric") %>%
  ggplot(., aes(x = PC1, y = PC2)) +
  geom_point(aes(color = LOCALITY, 
                 shape = paste(zea,symp_allo)), 
                 alpha = .8) + 
  labs(title = "PCA of mex.mex haplotypes at inv4m",
       x = paste("PC 1 -", round(e$values[1]), "%"),
       y = paste("PC 2 -", round(e$values[2]), "%"))
ggsave("../plots/PCA_inv4m_mex_haps.png", device = "png", 
       width = 12, height = 10, units = "in",
       dpi = 200)
# and just the maize haplotypes:
pca2 %>%
  filter(anc == "maize.maize" | hap_group == "maize_allopatric") %>%
  ggplot(., aes(x = PC1, y = PC2)) +
  geom_point(aes(color = LOCALITY, shape = paste(zea,symp_allo)), 
                 alpha = .8) + 
  labs(title = "PCA of maize.maize haplotypes at inv4m",
       x = paste("PC 1 -", round(e$values[1]), "%"),
       y = paste("PC 2 -", round(e$values[2]), "%"))
ggsave("../plots/PCA_inv4m_maize_haps.png", device = "png", 
       width = 12, height = 10, units = "in",
       dpi = 200)



# TO DO: write a function to make key plots for outlier regions
# also get coverage for this region, to plot too
# but for not plot pca for peak2 on chr4
m1 <- as.matrix(read.table("../data/outliers/chr4/peak2/PCA.cov", stringsAsFactors = F, header = F))
e1 <- eigen(m1)

# load inv4m_anc from plotLocalAncTracts.R for labels
pca0 = e1$vectors[, 1:2]
colnames(pca0) = paste0("PC", 1:2)
pca1 = cbind(pass1_allo4Low, pca0) %>%
  left_join(., peak2_anc, by = c("ID", "n"="hilo_id", "zea"))
# there are some NAs because they are from pops without ancestry calls
pca1$hap_group[is.na(pca1$hap_group)] <- paste(pca1$zea, pca1$symp_allo, sep = "_")[is.na(pca1$hap_group)]
pca1 %>%
  ggplot(., aes(x = PC1, y = PC2)) +
  geom_point(aes(color = hap_group))
pca1 %>%
  #filter(!hap_group == "maize_sympatric" & !hap_group == "mexicana_sympatric") %>%
  ggplot(., aes(x = PC1, y = PC2)) +
  geom_point(aes(color = anc, 
                 #size = top_anc_n,
                 size = est_coverage,
                 shape = paste(zea,symp_allo)),
             alpha = .8) + 
  labs(title = "PCA of haplotypes at peak 2 chr4",
       x = paste("PC 1 -", round(e1$values[1]), "%"),
       y = paste("PC 2 -", round(e1$values[2]), "%"))
ggsave("../plots/PCA_peak2.png", device = "png", 
       width = 12, height = 10, units = "in",
       dpi = 200)
# do the mexicana haplotypes appear to cluster at all by location?
# no, not really. I see this when I just plot mex.mex haplotypes
pca1 %>%
  filter(anc == "mex.mex" | hap_group == "mexicana_allopatric") %>%
  ggplot(., aes(x = PC1, y = PC2)) +
  geom_point(aes(color = LOCALITY, 
                 shape = paste(zea,symp_allo)), 
             alpha = .8) + 
  labs(title = "PCA of mex.mex haplotypes at peak2 chr4",
       x = paste("PC 1 -", round(e1$values[1]), "%"),
       y = paste("PC 2 -", round(e1$values[2]), "%"))
ggsave("../plots/PCA_peak2_mex_haps.png", device = "png", 
       width = 12, height = 10, units = "in",
       dpi = 200)
# and just the maize haplotypes:
pca1 %>%
  filter(anc == "maize.maize" | hap_group == "maize_allopatric") %>%
  ggplot(., aes(x = PC1, y = PC2)) +
  geom_point(aes(color = LOCALITY, shape = paste(zea,symp_allo)), 
             alpha = .8) + 
  labs(title = "PCA of maize.maize haplotypes at peak2 chr4",
       x = paste("PC 1 -", round(e1$values[1]), "%"),
       y = paste("PC 2 -", round(e1$values[2]), "%"))
ggsave("../plots/PCA_peak2_maize_haps.png", device = "png", 
       width = 12, height = 10, units = "in",
       dpi = 200)


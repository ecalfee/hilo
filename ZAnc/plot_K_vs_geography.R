#!/usr/bin/env Rscript
library(dplyr)
library(geosphere)
library(ade4)

# this script compares ancestry variance-covariance K matrix
# for maize or mexicana and geographic distance (from lat/lon)

# load variables from Snakefile
zea = snakemake@params[["zea"]]
# zea = "maize"
# zea = "mexicana"
K_file = snakemake@input[["K"]]
# K_file = paste0("ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/", zea, ".K.RData")
meta_file = snakemake@input[["meta_pop"]]
# meta_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/", zea, ".pop.meta.RData")


# load data
load(K_file)
load(meta_file)

# number of pops
n = nrow(K$K)
# calculate distance between each pair of pops
dist_matrix = matrix(NA, n, n)
colnames(dist_matrix) = colnames(K$K)
rownames(dist_matrix) = rownames(K$K)
for (i in 1:n){
  for (j in 1:n){
    dist_matrix[i, j] = distm(meta_pops[meta_pops$pop == colnames(dist_matrix)[i], 
                                        c("LON", "LAT")],
                              meta_pops[meta_pops$pop == rownames(dist_matrix)[j], 
                                        c("LON", "LAT")],
                              fun = distGeo)/1000 # distance in km    
  }
}

# use mantel test to calculate similarity between these two matrices
dist_anc = 1 - cov2cor(K$K)
#dist_anc = K$K
dist_anc[upper.tri(K$K, diag = F)] <- 0
dist_geo = dist_matrix
dist_geo[upper.tri(K$K, diag = T)] <- 0
ade4::mantel.rtest(as.dist(dist_anc), as.dist(dist_geo), nrepet = 999)
plot(dist_anc, dist_geo)
plot(K$K, dist_matrix)
cor(K$K[upper.tri(K$K, diag = F)], dist_matrix[upper.tri(K$K, diag = F)])
cor(dist_anc[lower.tri(K$K, diag = F)], dist_geo[lower.tri(K$K, diag = F)])

# how similar are maize and mexicana K matrices?



# plots PCA results
library(scales)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)# for color palette

# standard colors and labels for plots
blues = brewer.pal(n = 9, name = "Blues")[c(6,8)]
yellows = brewer.pal(n = 9, name = "YlOrBr")[c(4,6)]
colors_maize2mex = c(yellows, blues)
colors_alphabetical = colors_maize2mex[c(1,4,2,3)] # allo maize, allo mex, symp maize, symp mex


PREFIX <- "pass2_alloMAIZE"
#PREFIX <- "duplicates"
#PREFIX <- "hilo_alloMAIZE_MAIZE4LOW_RIMMA0625_small_and_duplicates"
PREFIX_metrics <- "hilo_alloMAIZE_MAIZE4LOW"
ids_dup = unlist(lapply(data.frame(ID = read.table(paste0("../samples/duplicates_IDs.list"), header = F, stringsAsFactors = F)$V1, stringsAsFactors = F),
                        function(x) rep(x, 2)))
ids_merged = read.table(paste0("../samples/", "hilo_alloMAIZE_MAIZE4LOW_RIMMA0625_small", "_IDs.list"), header = F, stringsAsFactors = F)$V1
#IDs = data.frame(ID = c(ids_merged, ids_dup))
#PREFIX <- "missing268-276_hilo_alloMAIZE_MAIZE4LOW"
IDs <- data.frame(ID = read.table(paste0("../samples/", PREFIX, "_IDs.list"), header = F, stringsAsFactors = F)$V1, stringsAsFactors = F)
metrics <- read.table(paste0("../filtered_bams/metrics/", PREFIX_metrics, ".flagstat.total"), header = F, stringsAsFactors = F)
colnames(metrics) <- c("ID", "total_reads_pass")
hilo <- read.table("../samples/hilo_meta.txt", stringsAsFactors = F, header = T, sep = "\t")
maize4low <- read.table("../samples/MAIZE4LOW_meta.txt", stringsAsFactors = F, header = T)
landraces <- read.table("../samples/alloMAIZE_meta.txt", stringsAsFactors = F, header = T)
seq_Jan2019 <- read.table("../data/HILO_raw_reads/Jan2019_IDs.list", stringsAsFactors = F, header = F)$V1

# elevation
pop_elev <- read.table("../data/riplasm/gps_and_elevation_for_sample_sites.txt",
                       stringsAsFactors = F, header = T, sep = "\t") %>%
  dplyr::select(., popN, ELEVATION)

# only for using RIMMA0625_small
#landraces$ID[landraces$ID == "RIMMA0625"] <- "RIMMA0625_small"

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
  left_join(., pop_elev, by = "popN") %>%
  mutate(., group = paste(symp_allo, zea, sep = "_"))
#meta$est_coverage[meta$ID=="RIMMA0625"] <- metrics$est_coverage[metrics$ID=="RIMMA0625_small"]
sum(meta$est_coverage < .025)

# which individuals are included in 'pass2'?
pass2 <- read.table("../samples/pass2_IDs.list", stringsAsFactors = F, header = F)$V1

# just plot numbers for pass2:
for (i in c(0, 0.25, 0.5, 1)) {
  p_counts <- meta %>%
    filter(., ID %in% pass2) %>%
    filter(., est_coverage >= i) %>%
    ggplot(aes(fill = group, x = LOCALITY)) +
    geom_histogram(stat = "count") +
    ggtitle(paste0("pass2; # individuals with coverage > ", i, "x")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggsave(paste0("plots/counts_pass2_morethan", i, "x_coverage.png"),
         plot = p_counts,
         device = "png",
         width = 12, height = 8, units = "in",
         dpi = 200)
}

# what is the total depth per population?
sum(meta$est_coverage) # total coverage
meta %>% # per pop/group
  filter(., ID %in% pass2) %>%
  group_by(., group, LOCALITY) %>%
  summarize(total_x = sum(est_coverage))
# plotted coverage per population:
meta %>%
  filter(., ID %in% pass2) %>%
  ggplot(aes(fill = group, x = reorder(LOCALITY, ELEVATION), y = est_coverage)) +
  geom_bar(stat = "identity") +
  ggtitle(paste0("pass2; total x coverage per group")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(paste0("plots/total_x_coverage_pass2.png"),
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)
# how evenly is that coverage distributed across individuals?
meta %>%
  filter(., ID %in% pass2) %>%
  ggplot(aes(color = group, x = reorder(LOCALITY, ELEVATION), y = est_coverage)) +
  geom_point() +
  ggtitle(paste0("pass2; individual x coverage")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(paste0("plots/ind_x_coverage_pass2_scatter_by_pop.png"),
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)


# how many individuals did we lose due to excluding HILO80 and PCR plate #2?
meta %>%
  filter(., group != "allopatric_maize") %>%
  filter(., est_coverage >= 0.25) %>%
  filter(., !(ID %in% pass2)) %>%
  ggplot(aes(fill = group, x = reorder(LOCALITY, ELEVATION))) +
  xlab("LOCALITY - (!) group/pop may not be accurate") +
  geom_histogram(stat = "count") +
  ggtitle("# individuals excluded with > 0.25x coverage") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("plots/counts_excluded_samples_with_morethan0.25x_coverage.png",
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)
meta %>%
  filter(., group != "allopatric_maize") %>%
  filter(., est_coverage >= 0.5) %>%
  filter(., !(ID %in% pass2)) %>%
  ggplot(aes(fill = group, x = reorder(LOCALITY, ELEVATION))) +
  xlab("LOCALITY - (!) group/pop may not be accurate") +
  geom_histogram(stat = "count") +
  ggtitle("# individuals excluded with > 0.5x coverage") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("plots/counts_excluded_samples_with_morethan0.5x_coverage.png",
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)

# what are the projections after resequencing individuals?
new_plants <- read.table("../samples/seeds_planted_5.1.19.csv", 
                         quote="",
                         sep = ",",
                         fill = F,
                         stringsAsFactors = F, header = T) %>%
  mutate(est_coverage = 1.5) %>% # projected coverage post-sequencing
  left_join(., unique(dplyr::select(meta, 
                                    c("popN", "zea", "symp_allo", "RI_ACCESSION", "GEOCTY", "LOCALITY", "ELEVATION", "group"))),
                      by = c("popN", "RI_ACCESSION"))
pass2_plus_reseq <- filter(meta, ID %in% pass2) %>%
  bind_rows(., new_plants) %>%
  mutate(., fam = substr(family, 1, 1)) %>%
  arrange(., -est_coverage) %>% # if duplicates, mark the lower coverage one as the duplicate
  mutate(., replaced = duplicated(dplyr::select(., c("popN", "fam")))) %>%
  mutate(., under0.5x = est_coverage < 0.5)

# about half of the individuals < 0.5x coverage are getting replaced
table(dplyr::select(pass2_plus_reseq, c("replaced", "under0.5x")))

# what's going on in plots?
pass2_plus_reseq %>%
  filter(LOCALITY=="Tenango del Aire") %>%
  group_by(zea) %>%
  summarize(., tot_cov = sum(est_coverage))

# what does projected coverage look like?
pass2_plus_reseq %>%
  filter(., !replaced) %>%
  filter(., (symp_allo == "allopatric" | !under0.5x)) %>%
  ggplot(aes(fill = group,  
             x = reorder(LOCALITY, ELEVATION), 
             y = est_coverage)) +
  geom_bar(stat = 'identity') +
  ggtitle(paste0("pass2 plus reseq; total x coverage per group")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(paste0("plots/total_x_coverage_pass2_plus_reseq.png"),
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)
# how evenly is that coverage distributed across individuals?
pass2_plus_reseq %>%
  filter(., !replaced) %>%
  ggplot(aes(color = group, x = reorder(LOCALITY, ELEVATION), y = est_coverage,
             shape = under0.5x)) +
  geom_point() +
  ggtitle(paste0("pass2 plus reseq; individual x coverage")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(paste0("plots/ind_x_coverage_pass2_plus_reseq_scatter_by_pop.png"),
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)
# how many individuals per population will be included in local ancestry inference?
pass2_plus_reseq %>%
  filter(., !replaced) %>%
  filter(., (symp_allo == "allopatric" | !under0.5x)) %>%
  ggplot(aes(fill = group,  
             x = reorder(LOCALITY, ELEVATION))) +
  geom_bar(stat = 'count') +
  ggtitle(paste0("pass2 plus reseq; samples >0.5x coverage per group")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(paste0("plots/counts_pass2_plus_reseq_morethan0.5x_coverage.png"),
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)

# table of population coverage
pass2_plus_reseq %>%
  filter(., !replaced) %>%
  filter(., (symp_allo == "allopatric" | !under0.5x)) %>%
  group_by(LOCALITY, zea) %>%
  summarise(., tot_cov = sum(est_coverage)) %>%
  summary(.)

# update: 32 plants didn't grow
failed <- read.table("../samples/2019_25_May_failed2grow.txt", stringsAsFactors = F, header = F, sep = " ")
colnames(failed) <- c("RIMMEA", "popN", "family")
pass2_plus_reseq$failed <- !pass2_plus_reseq$replaced & paste0(pass2_plus_reseq$popN, substr(pass2_plus_reseq$family, 1, 1)) %in% paste0(failed$popN, substr(failed$family, 1, 1))
# now how many per pop:
pass2_plus_reseq %>%
  filter(., !replaced) %>%
  filter(., (symp_allo == "allopatric" | !under0.5x)) %>%
  ggplot(aes(fill = failed,  
             x = reorder(LOCALITY, ELEVATION))) +
  geom_bar(stat = 'count') +
  facet_wrap(~group) +
  ggtitle(paste0("pass2 plus reseq; samples >0.5x coverage per group")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(paste0("plots/counts_pass2_plus_reseq_morethan0.5x_coverage_failed_grow2.png"),
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)
# coverage
pass2_plus_reseq %>%
  filter(., !replaced) %>%
  filter(., (symp_allo == "allopatric" | !under0.5x)) %>%
  ggplot(aes(fill = failed,  
             x = reorder(LOCALITY, ELEVATION), 
             y = est_coverage)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~group) +
  ggtitle(paste0("pass2 plus reseq; total x coverage per group")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(paste0("plots/total_x_coverage_pass2_plus_reseq_failed_grow2.png"),
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)





# get PCA data
cov_dir <- file.path("results/PCA", PREFIX)
cov_data <- read.table(file.path(cov_dir, "whole_genome.cov"),
                       header = F, stringsAsFactors = F)
# PC's are each column of dataframe pca:
# i.e. PC1 is V1, PC2 is V2 etc.
pca <- eigen(cov_data) # take PCA of covariance matrix
n = 4 # make small dataframe w/ only first n eigenvectors
pca_small <- data.frame(pca$vectors[ , 1:n])
colnames(pca_small) = paste0("PC", 1:n)
# rounded eigen values
PC_var_explained = round(pca$values, 2)

#IDs <- data.frame(ID = sapply(IDs, function(x) rep(x, 2)))
# quick plot of new data:
d <- bind_cols(IDs, pca_small) %>%
  left_join(., meta, by = "ID")
d %>%
  filter(., ID %in% seq_Jan2019 | (symp_allo == "allopatric" & zea == "maize")) %>%
  ggplot(., aes(x = -PC1, y = PC2, color = group,
                size = est_coverage)) +
  geom_point(alpha = .5) +
  scale_colour_manual(values = colors_alphabetical) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  ggtitle("newly seq. HILO Jan19 all; SNPs spaced >.01cM")

ggsave("plots/PC_new_HILO_seq_Jan19.png",
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)

d %>%
  filter(., (ID %in% seq_Jan2019 & n > 200) | (symp_allo == "allopatric" & zea == "maize")) %>%
  ggplot(., aes(x = -PC1, y = PC2, color = group,
                size = est_coverage)) +
  geom_point(alpha = .5) +
  scale_colour_manual(values = colors_alphabetical[c(1,3:4)]) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  ggtitle("newly seq. HILO Jan19 only (n>200); SNPs spaced >.01cM")
ggsave("plots/PC_new_HILO_seq_post_Jan19_only.png",
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)


d %>%
  filter(., !(ID %in% seq_Jan2019) | (symp_allo == "allopatric" & zea == "maize")) %>%
  ggplot(., aes(x = -PC1, y = PC2, color = group,
                size = est_coverage)) +
  geom_point(alpha = .75) +
  scale_colour_manual(values = colors_alphabetical[c(1,3:4)]) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  ggtitle("old seq. HILO only; SNPs spaced >.01cM")
ggsave("plots/PC_old_HILO_seq_before_Jan19_only.png",
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)

d %>%
  filter(., (ID %in% seq_Jan2019 & n <= 200) | (symp_allo == "allopatric" & zea == "maize")) %>%
  #filter(., (ID %in% seq_Jan2019 & n <= 200)) %>%
  ggplot(., aes(x = -PC1, y = PC2, color = group,
                size = est_coverage)) +
  scale_colour_manual(values = colors_alphabetical) +
  geom_point(alpha = .75) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC1 (", PC_var_explained[2], "%)")) +
  ggtitle("new seq. HILO, overlap w/ pass1 only; SNPs spaced >.01cM")
ggsave("plots/PC_HILO_seq_Jan19_overlap_w_pass1.png",
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)

# just look at allopatric maize
d %>%
  #filter(., (symp_allo == "allopatric" & zea == "maize")) %>%
  ggplot(., aes(x = -PC1, y = PC2,
                color = ifelse((symp_allo == "allopatric" & zea == "maize"),
                               LOCALITY, "hilo_samples"),
                size = est_coverage)) +
  geom_point(alpha = .5) +
  #scale_colour_manual(values = colors_alphabetical) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  ggtitle("structure within allopatric maize; SNPs spaced >.01cM")
ggsave("plots/PC_allopatric_maize_landraces.png",
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)
# allopatric maize PC3 & 4
d %>%
  ggplot(., aes(x = PC3, y = PC4,
                color = ifelse((symp_allo == "allopatric" & zea == "maize"),
                               LOCALITY, "hilo_samples"),
                size = est_coverage)) +
  geom_point(alpha = .5) +
  #scale_colour_manual(values = colors_alphabetical) +
  xlab(paste0("PC3 (", PC_var_explained[3], "%)")) +
  ylab(paste0("PC4 (", PC_var_explained[4], "%)")) +
  ggtitle("PC3 & 4 within allopatric maize; SNPs spaced >.01cM")
ggsave("plots/PC34_allopatric_maize_landraces.png",
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)

# exclude Anne's allopatric 4 lowland populations:
pca2 <- eigen(cov_data[1:291,1:291]) # take PCA of covariance matrix
n = 4 # make small dataframe w/ only first n eigenvectors
pca_small2 <- data.frame(pca2$vectors[ , 1:n])
colnames(pca_small2) = paste0("PC", 1:n)
# rounded eigen values
PC_var_explained2 = round(pca2$values, 2)
#1:291
d2 <- data.frame(ID = IDs[1:291,], pca_small2, stringsAsFactors = F) %>%
  left_join(., meta, by = "ID")
d2 %>%
  ggplot(., aes(x = -PC1, y = PC2, color = group,
                size = est_coverage)) +
  geom_point(alpha = .5) +
  #scale_colour_manual(values = colors_alphabetical) +
  xlab(paste0("PC1 (", PC_var_explained2[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained2[2], "%)")) +
  ggtitle("Li's landraces and all hilo samples; SNPs spaced >.01cM")
ggsave("plots/PCA_Lis_landraces_plus_hilo.png",
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)



# old data (before 2019 sequencing)
cov_data_old <- read.table("../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedPCA/PCA/whole_genome.cov",
                           header = F, stringsAsFactors = F)

# PC's are each column of dataframe pca:
# i.e. PC1 is V1, PC2 is V2 etc.
pca_old <- eigen(cov_data_old) # take PCA of covariance matrix
pca_small_old <- data.frame(pca_old$vectors[ , 1:n])
colnames(pca_small_old) = paste0("PC", 1:n)
# rounded eigen values
PC_var_explained_old = round(pca_old$values, 2)
pass1_allo4Low <- read.table("../data/pass1_allo4Low_ids.txt", stringsAsFactors = F,
                             header = T, sep = "\t")
# join bams and first 10 PCs of covariance PCA data by position (CAUTION - bam list order and admix results MUST MATCH!)
d_old <- bind_cols(pass1_allo4Low, pca_small_old)  %>%
  #arrange(., popN) %>%
  #arrange(., zea) %>%
  #arrange(., symp_allo) %>%
  mutate(., group = paste(symp_allo, zea, sep = "_"))

# plot first PC's
# PC1 and 2, all no filter
p12_coverage = d_old %>%
  ggplot(., aes(PC1, PC2, alpha = .5)) +
  xlab(paste0("PC1 (", PC_var_explained_old[1], "%)")) +
  ylab(paste0("PC1 (", PC_var_explained_old[2], "%)")) +
  ggtitle("PC1 and PC2 (all ind's); SNPs spaced >.01cM") +
  geom_point(aes(color = group, size = est_coverage, shape = (ID == "HILO80"))) +
  scale_colour_manual(values = colors_alphabetical)
plot(p12_coverage)
ggsave("plots/PCA_12_by_coverage.png", plot = p12_coverage, device = png(),
       width = 12, height = 8, units = "in",
       dpi = 200)
# PC1 and 2, excluding super low coverae ind's
p12 = d_old %>%
  filter(., est_coverage >= .05) %>%
  ggplot(., aes(PC1, PC2)) +
  ggtitle("PC1 and PC2 (excl. low coverage ind's < .05x); SNPs spaced >.01cM") +
  xlab(paste0("PC1 (", PC_var_explained_old[1], "%)")) +
  ylab(paste0("PC1 (", PC_var_explained_old[2], "%)")) +
  geom_point(aes(color = group)) +
  scale_colour_manual(values = colors_alphabetical)
plot(p12)
ggsave("plots/PCA_12.png", plot = p12, device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)
# PC3 and PC4, excluding low coverage ind's
p34 = d_old %>%
  filter(., est_coverage >= .05) %>%
  ggplot(., aes(PC3, PC4)) +
  ggtitle("PC3 and PC4 (excl. low coverage ind's < .05x); SNPs spaced >.01cM") +
  xlab(paste0("PC3 (", PC_var_explained_old[3], "%)")) +
  ylab(paste0("PC4 (", PC_var_explained_old[4], "%)")) +
  geom_point(aes(color = group)) +
  scale_colour_manual(values = colors_alphabetical)
plot(p34)
ggsave("plots/PCA_34.png", plot = p34, device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)

# Use vars() to supply faceting variables:
p12 + facet_wrap(~LOCALITY)
ggsave("plots/PCA_facet_by_population.png",
       plot = p12 + facet_wrap(~LOCALITY), device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)


# repeat same analyses for SNP data pruned every 1000 SNPs
cov_data_1000 <- read.table("../data/geno_lik/merged_pass1_all_alloMaize4Low_16/prunedBy1000/PCA/whole_genome.cov",
                           header = F, stringsAsFactors = F)

pca_1000 <- eigen(cov_data_1000) # take PCA of covariance matrix
n = 4 # make small dataframe w/ only first n eigenvectors
pca_small_1000 <- data.frame(pca_1000$vectors[ , 1:n])
colnames(pca_small_1000) = paste0("PC", 1:n)

# join bams and first 10 PCs of covariance PCA data by position (CAUTION - bam list order and admix results MUST MATCH!)
d_1000 <- bind_cols(pass1_allo4Low, pca_small_1000)  %>%
  #arrange(., popN) %>%
  #arrange(., zea) %>%
  #arrange(., symp_allo) %>%
  mutate(., group = paste(symp_allo, zea, sep = "_"))

# rounded eigen values
PC_var_explained_1000 = round(pca_1000$values, 2)

# plot first PC's
# PC1 and 2, all no filter
p12_coverage_1000 = d_1000 %>%
  ggplot(., aes(PC1, PC2, alpha = .5)) +
  xlab(paste0("PC1 (", PC_var_explained_1000[1], "%)")) +
  ylab(paste0("PC1 (", PC_var_explained_1000[2], "%)")) +
  ggtitle("PC1 and PC2 (all ind's); every 1000th SNP") +
  geom_point(aes(color = group, size = est_coverage)) +
  scale_colour_manual(values = colors_alphabetical)
ggsave("plots/PCA_12_by_coverage_pruned1000.png",
       plot = p12_coverage_1000, device = png(),
       width = 12, height = 8, units = "in",
       dpi = 200)
# PC1 and 2, excluding super low coverae ind's
p12_1000 = d_1000 %>%
  filter(., est_coverage >= .05) %>%
  ggplot(., aes(PC1, PC2)) +
  ggtitle("PC1 and PC2 (excl. low coverage ind's < .05x; every 1000th SNP)") +
  xlab(paste0("PC1 (", PC_var_explained_1000[1], "%)")) +
  ylab(paste0("PC1 (", PC_var_explained_1000[2], "%)")) +
  geom_point(aes(color = group)) +
  scale_colour_manual(values = colors_alphabetical)
plot(p12_1000)
ggsave("plots/PCA_12_pruned1000.png", plot = p12_1000, device = png(),
       width = 12, height = 8, units = "in",
       dpi = 200)
# PC3 and PC4, excluding low coverage ind's
p34_1000 = d_1000 %>%
  filter(., est_coverage >= .05) %>%
  ggplot(., aes(PC3, PC4)) +
  ggtitle("PC3 and PC4 (excl. low coverage ind's < .05x); every 1000th SNP") +
  xlab(paste0("PC3 (", PC_var_explained_1000[3], "%)")) +
  ylab(paste0("PC4 (", PC_var_explained_1000[4], "%)")) +
  geom_point(aes(color = group)) +
  scale_colour_manual(values = colors_alphabetical)
plot(p34_1000)
ggsave("plots/PCA_34_pruned1000.png", plot = p34_1000, device = png(),
       width = 12, height = 8, units = "in",
       dpi = 200)

# Use vars() to supply faceting variables:
p12_1000 + facet_wrap(~LOCALITY)
ggsave("plots/PCA_facet_by_population_pruned1000.png",
       plot = p12_1000 + facet_wrap(~LOCALITY), device = png(),
       width = 12, height = 8, units = "in",
       dpi = 200)


# trouble shooting new vs old sequencing:
# re-plot PCA, but by facet_wrap(~lane)

d %>%
  left_join(., jan19, by = "ID") %>%
  mutate(., lane2 = ifelse(is.na(novogene_ID), "pass1_only", lane)) %>%
  ggplot(., aes(PC1, PC2, alpha = .5)) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  ggtitle("new PCA facet by seq lane") +
  geom_point(aes(color = group, size = est_coverage)) +
  scale_colour_manual(values = colors_alphabetical) +
  facet_wrap(~lane2)
ggsave("plots/jan2019_new_PCA_facet_by_seq_lane.png",
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)



# compare PCA positions for individuals resequenced in Jan 2019
par(mfrow=c(1,2))
# old data
d_1000 %>%
  filter(., ID %in% seq_Jan2019 & ID %in% pass1_allo4Low$ID) %>%
  ggplot(., aes(-PC1, PC2, alpha = .5)) +
  xlab(paste0("PC1 (", PC_var_explained_1000[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained_1000[2], "%)")) +
  ggtitle("PC1 and PC2 - old PCA (ind's that overlap with seq2019)") +
  geom_point(aes(color = ID, size = est_coverage)) +
  facet_wrap(~group)
ggsave("plots/jan2019_old_PCA_overlap.png",
       device = "png",
       width = 20, height = 6, units = "in",
       dpi = 200)

# new data
d %>%
  filter(ID %in% seq_Jan2019 & ID %in% pass1_allo4Low$ID) %>%
  ggplot(., aes(PC1, PC2, alpha = .5)) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  ggtitle("PC1 & PC2 - new seq 2019 (ind's that overlap with pass1)") +
  geom_point(aes(color = ID, size = est_coverage)) +
  facet_wrap(~group)
ggsave("plots/jan2019_new_PCA_overlap.png",
       device = "png",
       width = 20, height = 6, units = "in",
       dpi = 200)

# testing with new population labels
l678_raw <- read.table("../samples/HILO_test_lane6_7_8.csv",
                 sep = ",", header = F, stringsAsFactors = F)
colnames(l678_raw) <- c("ID", "popN_fam")
l678 <- l678_raw %>%
  separate(data = ., col = popN_fam, sep = "_", c("popN", "family"), extra = "merge") %>%
  mutate(popN = as.integer(popN))

# what population assignments changed? Basically HILO202 - HILO217 get reassigned
left_join(l678, d, by = "ID") %>%
  filter(., popN.x!=popN.y) %>%
  select(ID, popN.x, popN.y)
# get metadata for updated population labels
d678 <- d %>%
  select(ID, PC1, PC2, n) %>%
  left_join(., l678, by = "ID") %>%
  left_join(., unique(meta[ , c("popN", "zea", "symp_allo", "RI_ACCESSION", "GEOCTY", "LOCALITY", "group")]), by = "popN")
# plot with updated population/group labels
d678 %>%
  left_join(., select(jan19, c("ID", "lane")), by = "ID") %>%
  filter(!is.na(popN)) %>%
  ggplot(., aes(PC1, PC2, alpha = .5)) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  ggtitle("PC1 & PC2 - pools 6-8, new family ID's)") +
  geom_point(aes(color = group, shape = (n <= 219))) +
  facet_wrap(~lane)
ggsave("plots/jan2019_new_family_IDs_PCA.png",
       device = "png",
       width = 20, height = 6, units = "in",
       dpi = 200)

d %>%
  left_join(., select(jan19, c("ID", "lane")), by = "ID") %>%
  filter(!is.na(popN)) %>%
  ggplot(., aes(PC1, PC2, alpha = .5)) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  ggtitle("new PCA family Ids facet by seq lane") +
  geom_point(aes(color = group, size = est_coverage)) +
  scale_colour_manual(values = colors_alphabetical) +
  facet_wrap(~lane)
ggsave("plots/jan2019_new_family_Ids_PCA_facet_by_seq_lane.png",
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)

# are there individuals with high enough coverage in both seq runs
# to compare they are the same individual?
metrics_old <- read.table(paste0("../data/filtered_bam/pass1.all.metrics.raw.Nreads"), header = T, stringsAsFactors = F)

# get lane data
jan19_lanes <- read.table("../data/HILO_raw_reads/Jan2019_lanes.list", stringsAsFactors = F)$V1
jan19_IDs <- read.table("../data/HILO_raw_reads/Jan2019_IDs.list", stringsAsFactors = F)$V1
jan19 <- read.table("../data/HILO_raw_reads/Jan2019_id_link_HILO_novogene.txt", stringsAsFactors = F)
colnames(jan19) <- c("ID", "novogene_ID", "Pop_Fam")
#table(jan19$ID == jan19_IDs) # matches, great, same order
jan19$lane <- jan19_lanes
#table(jan19$lane, substr(jan19$novogene_ID, 1, 5)) # also matches, good

# read in adapter and plate data
pools <- read.table("../samples/hilo_plates.txt", header = T, stringsAsFactors = F)



# plot updated PCA now with diff. shapes for diff. pcr plates
d_updated <- bind_cols(IDs, pca_small) %>%
  left_join(., metrics, by = "ID") %>%
  left_join(., dplyr::select(pools, c("ID", "popN", "lane", "plate_number")), by = "ID") %>%
  left_join(., unique(dplyr::select(meta, c("popN", "zea", "symp_allo", "RI_ACCESSION", "GEOCTY", "LOCALITY", "group"))), by = "popN") %>%
  filter(., group != "allopatric_maize") %>% # do separately
  mutate(., n = as.integer(substr(ID, 5, 8))) %>%
  bind_rows(., filter(d, group == "allopatric_maize")) %>% # add in metadata for allopatric maize
  mutate(., plate_number_pcr = ifelse(group == "allopatric_maize", "external_seq", plate_number))

d_updated %>%
  ggplot(., aes(PC1, PC2, alpha = .5)) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  ggtitle("new PCA - show pcr plates - facet by seq lane") +
  geom_point(aes(color = group, size = est_coverage, shape = plate_number_pcr)) +
  scale_colour_manual(values = colors_alphabetical) +
  facet_wrap(~lane)
ggsave("plots/jan2019_new_seq_by_pcr_plate_and_novogene_lane.png",
       device = "png",
       width = 8, height = 6, units = "in",
       dpi = 200)

d_updated %>%
  #filter(., ID %in% seq_Jan2019 & ID %in% pass1_allo4Low$ID | group == "allopatric_maize") %>%
  filter(., duplicated(ID) | group == "allopatric_maize") %>%
  ggplot(., aes(PC1, PC2, alpha = .5)) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  ggtitle("PCA of samples with repeated sequencing") +
  geom_point(aes(color = group, size = est_coverage, shape = plate_number_pcr)) +
  scale_colour_manual(values = colors_alphabetical) +
  facet_wrap(~lane)
ggsave("plots/jan2019_are_duplicates_new_seq_by_pcr_plate_and_novogene_lane.png",
       device = "png",
       width = 8, height = 6, units = "in",
       dpi = 200)

# how many samples look like they're on the wrong side of the PCA for maize/mexicana ancestry?
d_updated %>%
  filter(., ((zea == "maize" & PC1 < 0) | (zea == "mexicana" & PC1 > 0)) & plate_number_pcr != "plate2") %>%
  View(.) # just a couple outside of plate 2
# plot these samples
d_updated %>%
  filter(., ((zea == "maize" & PC1 < 0) | (zea == "mexicana" & PC1 > 0))) %>%
  ggplot(., aes(PC1, PC2)) +
  xlab(paste0("mexicana like to maize like - PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  ggtitle("samples where PCA positions don't match mex/maize label") +
  geom_point(aes(color = group, size = est_coverage), alpha = .5) +
  scale_colour_manual(values = colors_alphabetical[3:4]) +
  facet_wrap(~plate_number_pcr) +
  geom_text(aes(label=ifelse(plate_number_pcr != "plate2",
                             ID, '')), hjust=0, vjust=0)
ggsave("plots/jan2019_samples_wrong_side_of_PC1.png",
       device = "png",
       width = 8, height = 6, units = "in",
       dpi = 200)
# looks at data frame for these samples:
d_updated %>%
  filter(., ((zea == "maize" & PC1 < 0) | (zea == "mexicana" & PC1 > 0))) %>%
  dplyr::select(ID, PC1, PC2, est_coverage, lane, plate_number, group, LOCALITY, popN) %>%
  mutate(label_switch = ifelse(abs(PC1) > 0.01, "likely", "possible")) %>%
  write.table(., "../samples/possible_label_switches_2019_hilo.txt",
              quote = F, col.names = T, row.names = F, sep = "\t")


# coverage looks particularly poor for lane 7;
# a lot higher % of libraries failed for lanes 6 & 7:
table(d_updated$lane, d_updated$est_coverage <= .05)


# rerun PCAngsd with separate entries for bams of the same individual,
# but before merging across sequencing pools to confirm it's the same sample
# plot duplicates:
# ok I don't understand why HILO16 is the only one driving the PCA at all
d %>%
  filter(ID != "HILO16") %>%
  ggplot(., aes(x = PC1, y = PC2, color = ID)) +
  geom_point(alpha = .5) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  ggtitle("duplicates") #+
  #geom_text(aes(label=ID), hjust=0, vjust=0)

# this is pointless - I need to include a broader range of samples in my PCA and then see what's happening with these.
pca_no16 <- eigen(cov_data[IDs$ID != "HILO16", IDs$ID != "HILO16"])
#identify(pca_no16$vectors[,1], pca_no16$vectors[,2])

# big picture with duplicates
# combine GL files to look at duplicates of the same samples along with their merged files on a PCA:
d1 <- read.table(gzfile("results/thinnedSNPs/hilo_alloMAIZE_MAIZE4LOW_RIMMA0625_small/whole_genome.beagle.gz"), header = T, stringsAsFactors = F)
d2 <- read.table(gzfile("results/thinnedSNPs/duplicates/whole_genome.beagle.gz"), header = T, stringsAsFactors = F)
# I'm ignoring SNPs where my duplicates had no data (i.e. not in d2)
d3 <- right_join(d1, d2, by = c("marker", "allele1", "allele2"))
colnames(d3) <- d3_header
d3_header <- c("marker", "allele1", "allele2", paste0("Ind", sapply(0:(dim(d3)[2]/3-2), function(x) rep(x, 3))))
d3_IDs <- data.frame(ID = c(read.table(paste0("../samples/", "hilo_alloMAIZE_MAIZE4LOW_RIMMA0625_small", "_IDs.list"), header = F, stringsAsFactors = F)$V1,
                            sapply(read.table(paste0("../samples/", "duplicates", "_IDs.list"), header = F, stringsAsFactors = F)$V1,
                                   function(x) rep(x, 2))),
                     source = c(rep("merged", 306), rep(c("dup1", "dup2"), 37)),
                     bam = c(read.table(paste0("../samples/", "hilo_alloMAIZE_MAIZE4LOW_RIMMA0625_small", "_bams.list"), header = F, stringsAsFactors = F)$V1,
                             read.table(paste0("../samples/", "duplicates", "_bams.list"), header = F, stringsAsFactors = F)$V1))
gz1 <- gzfile("results/thinnedSNPs/hilo_alloMAIZE_MAIZE4LOW_RIMMA0625_small_and_duplicates/whole_genome.beagle.gz", "w")
write.table(as.data.frame(d3, col.names = d3_header, check.names = F), # slow
            gz1,
            sep = "\t", quote = F,
            row.names = F, col.names = T)
close(gz1)

# need to fix merge:
d %>%
  ggplot(., aes(x = -PC1, y = PC2, color = group, size = est_coverage, shape = (ID %in% ids_dup))) +
  geom_point(alpha = .5) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  ggtitle("duplicates, merged data, and hilo")
ggsave("plots/resequenced_duplicates_PCA_perspective_with_other_hilo.png",
       device = "png",
       width = 8, height = 8, units = "in",
       dpi = 200)

# zoomed in on just the duplicates
d %>%
  filter(ID %in% ids_dup) %>%
  filter(duplicated(ID)) %>% # only look at original (unmerged)
  ggplot(., aes(x = -PC1, y = PC2, color = ID, shape = group, size = est_coverage)) +
  geom_point(alpha = .5) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  ggtitle("duplicates")
ggsave("plots/resequenced_duplicates_PCA.png",
       device = "png",
       width = 8, height = 8, units = "in",
       dpi = 200)
d %>%
  filter(ID %in% ids_dup) %>%
  filter(duplicated(ID)) %>% # only look at original (unmerged)
  ggplot(., aes(x = -PC1, y = PC2, color = ID, shape = group)) +
  geom_point(alpha = .5) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  ggtitle("duplicates") +
  facet_wrap(~ID)
ggsave("plots/resequenced_duplicates_PCA_facet_ID.png",
       device = "png",
       width = 15, height = 15, units = "in",
       dpi = 200)


# plot PCA or pass2_alloMAIZE individuals:
d %>%
  mutate(depth = ifelse(group == "allopatric_maize", 2, est_coverage)) %>%
  ggplot(., aes(x = PC1, y = PC2)) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  geom_point(alpha = 0.5, aes(color = group, size = depth)) +
  scale_colour_manual(values = colors_alphabetical) +
  ggtitle(paste("PCA HiLo and reference ind's colored by maize/mexicana group"))
ggsave(paste0("plots/PCA_", PREFIX, ".png"),
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)
d %>%
  mutate(depth = ifelse(group == "allopatric_maize", 2, est_coverage)) %>%
  ggplot(., aes(x = PC1, y = PC2)) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  geom_point(alpha = 0.5, aes(color = LOCALITY, size = depth)) +
  ggtitle(paste("PCA HiLo and reference ind's colored by population"))
ggsave(paste0("plots/PCA_", PREFIX, "_colored_by_pop.png"),
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)

d %>%
  mutate(depth = ifelse(group == "allopatric_maize", 2, est_coverage)) %>%
  ggplot(., aes(x = PC1, y = PC2)) +
  xlab(paste0("PC1 (", PC_var_explained[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained[2], "%)")) +
  geom_point(aes(color = ELEVATION, size = depth)) + # coverage truncated to be viewable
  ggtitle(paste("PCA HiLo colored by elevation")) +
  scale_color_gradientn(colors = brewer.pal(7, "YlGnBu")) # change the colors
ggsave(paste0("plots/PCA_", PREFIX, "_colored_by_elevation.png"),
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)

# what does the PCA look like for just inv4m, a large inversion on chr4?
# get PCA data
cov_dir_inv4m <- paste0("results/PCA/", PREFIX, "_inv4m")
cov_data_inv4m <- read.table(file.path(cov_dir_inv4m, "whole_genome.cov"),
                       header = F, stringsAsFactors = F)
# PC's are each column of dataframe pca:
# i.e. PC1 is V1, PC2 is V2 etc.
pca_inv4m <- eigen(cov_data_inv4m) # take PCA of covariance matrix
n = 4 # make small dataframe w/ only first n eigenvectors
pca_small_inv4m <- data.frame(pca_inv4m$vectors[ , 1:n])
colnames(pca_small_inv4m) = paste0("PC", 1:n)
# rounded eigen values
PC_var_explained_inv4m = round(pca_inv4m$values, 2)

#IDs <- data.frame(ID = sapply(IDs, function(x) rep(x, 2)))
# quick plot of new data:
d_inv4m <- bind_cols(IDs, pca_small_inv4m) %>%
  left_join(., meta, by = "ID")
d_inv4m %>%
  mutate(depth = ifelse(group == "allopatric_maize", 2, est_coverage)) %>%
  ggplot(., aes(x = PC1, y = PC2)) +
  xlab(paste0("PC1 (", PC_var_explained_inv4m[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained_inv4m[2], "%)")) +
  geom_point(aes(color = LOCALITY, shape = group, size = depth)) +
  ggtitle(paste("PCA HiLo inv4m locus")) #+
  #scale_color_gradientn(colors = brewer.pal(7, "YlGnBu")) # change the colors
ggsave(paste0("plots/PCA_", PREFIX, "inv4m.png"),
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)
d_inv4m %>%
  mutate(depth = ifelse(group == "allopatric_maize", 2, est_coverage)) %>%
  ggplot(., aes(x = PC1, y = PC2)) +
  xlab(paste0("PC1 (", PC_var_explained_inv4m[1], "%)")) +
  ylab(paste0("PC2 (", PC_var_explained_inv4m[2], "%)")) +
  geom_point(aes(color = ELEVATION, shape = group, size = depth)) +
  ggtitle(paste("PCA HiLo inv4m locus")) +
scale_color_gradientn(colors = brewer.pal(7, "YlGnBu")) # change the colors
ggsave(paste0("plots/PCA_", PREFIX, "inv4m_colored_by_elevation.png"),
       device = "png",
       width = 12, height = 8, units = "in",
       dpi = 200)
d_inv4m %>%
  mutate(depth = ifelse(group == "allopatric_maize", 2, est_coverage)) %>%
  ggplot(., aes(x = PC3, y = PC4)) +
  xlab(paste0("PC3 (", PC_var_explained_inv4m[3], "%)")) +
  ylab(paste0("PC4 (", PC_var_explained_inv4m[4], "%)")) +
  geom_point(aes(color = LOCALITY, shape = group, size = depth)) +
  ggtitle(paste("PCA HiLo inv4m locus - PC 3 & 4"))



# TO DO: possibly rerun PCAngsd with separate entries for bams of the same individual,
# but before merging across sequencing pools to confirm it's the same sample
# (include all final merged individuals too -- just add the pre-merging individual copies too)

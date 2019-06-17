# In this script I look at ancestry around genes associated with pathways that may be under strong selection vs. other genes and the background ancestry
# to test hypotheses:
# (I) Domestication genes are enriched in regions with low introgression
# (II) Barrier candidate genes are segregating in these pops. Maybe fewer barrier genes in pops with more gene flow. (predictions aren't super clear)
# (III) Mexicana ancestry has contributed to genetic variation selected on for flowering time. So I expect an enrichment of flowering time genes in positive
# zAnc ~ zElev hist (higher mexicana freq. at higher elevations)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)


# plot ancestry at candidate incompatibility loci:
GRMZM2G410783 <- read.table("../data/incompatibility_loci/GRMZM2G410783_v4_coord.bed", stringsAsFactors = F, sep = "\t")
d$GRMZM2G410783 <- (d$chr == GRMZM2G410783$V1 & d$pos >= GRMZM2G410783$V2 & d$pos <= GRMZM2G410783$V3)
AC231426.1_FG002 <- read.table("../data/incompatibility_loci/AC231426.1_FG002_v4_coord.bed", stringsAsFactors = F, sep = "\t")
d$AC231426.1_FG002 <- (d$chr == AC231426.1_FG002$V1 & d$pos >= AC231426.1_FG002$V2 & d$pos <= AC231426.1_FG002$V3)
d %>%
  ggplot(aes(x = GRMZM2G410783, y = pop_meanAlpha_mex)) +
  geom_boxplot()
d %>%
  ggplot(aes(x = AC231426.1_FG002, y = pop_meanAlpha_mex)) +
  geom_boxplot()
# local ancestry around loci of interest
d %>%
  tidyr::gather(., "zea", "meanMexAnc", c(pop_meanAlpha_maize, pop_meanAlpha_mex)) %>%
  filter(chr == 5 & pos >= 150000000 & pos <= 155000000) %>%
  tidyr::gather(., "locus", "within_locus", c(GRMZM2G410783, AC231426.1_FG002)) %>%
  #ggplot(aes(x = pos, y = meanMexAnc, color = GRMZM2G410783 | AC231426.1_FG002, shape = zea)) +
  ggplot(aes(x = pos, y = meanMexAnc, color = within_locus, shape = zea)) +
  geom_point() +
  facet_wrap(~locus) +
  ggtitle("mean mexicana ancestry in maize and mexicana at putative Ga2 loci")
ggsave("plots/mex_freq_around_candidate_incompatibility_loci_combined_pops.png",
       height = 6, width = 8, units = "in", device = "png")
# zoomed out view with boxplot
d %>%
  rename(maize = pop_meanAlpha_maize) %>%
  rename(mexicana = pop_meanAlpha_mex) %>%
  tidyr::gather(., "zea", "meanMexAnc", c(maize, mexicana)) %>%
  tidyr::gather(., "locus", "within_locus", c(GRMZM2G410783, AC231426.1_FG002)) %>%
  #mutate(zea = sapply(zea, function(x) strsplit(x, split = "_")[[1]][3])) %>%
  ggplot(aes(x = zea, y = meanMexAnc, color = within_locus)) +
  geom_boxplot() +
  facet_wrap(~locus) +
  ggtitle("mean mexicana ancestry in maize and mexicana at putative Ga2 loci compared to genomewide") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/mex_freq_boxplot_genomewide_vs_candidate_incompatibility_loci_combined_pops.png",
       height = 6, width = 8, units = "in", device = "png")

# plot by population:
# at female locus
bind_cols(d, all_anc) %>%
  tidyr::gather(., "pop", "meanMexAnc", meta.pops$pop) %>%
  left_join(., meta.pops, by = "pop") %>%
  dplyr::group_by(., pop, zea, LOCALITY, ELEVATION, GRMZM2G410783) %>%
  summarise(meanMexAnc = mean(meanMexAnc)) %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = meanMexAnc, color = GRMZM2G410783)) +
  geom_point() +
  facet_wrap(~zea) +
  ggtitle("mean mexicana ancestry in maize and mexicana at GRMZM2G410783 compared to genomewide") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/mex_freq_boxplot_genomewide_vs_candidate_incompatibility_GRMZM2G410783_ind_pops.png",
       height = 6, width = 8, units = "in", device = "png")
# male locus:
bind_cols(d, all_anc) %>%
  tidyr::gather(., "pop", "meanMexAnc", meta.pops$pop) %>%
  left_join(., meta.pops, by = "pop") %>%
  dplyr::group_by(., pop, zea, LOCALITY, ELEVATION, AC231426.1_FG002) %>%
  summarise(meanMexAnc = mean(meanMexAnc)) %>%
  ggplot(aes(x = reorder(LOCALITY, ELEVATION), y = meanMexAnc, color = AC231426.1_FG002)) +
  geom_point() +
  facet_wrap(~zea) +
  ggtitle("mean mexicana ancestry in maize and mexicana at AC231426.1_FG002 compared to genomewide") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/mex_freq_boxplot_genomewide_vs_candidate_incompatibility_AC231426.1_FG002_ind_pops.png",
       height = 6, width = 8, units = "in", device = "png")

# mexicana ancestry around domestication genes
domestication_genes <- read.table("../data/domestication/gene_model_translation_to_APGv4_Zm00001d.2_hits.txt", stringsAsFactors = F)$V1
domestication <- data.frame(ID = read.table("../data/domestication/Hufford_2012_domestication_genes_names_only.csv", stringsAsFactors = F, header = F)$V1,
                            stringsAsFactors = F)
flowering_pathway <- read.table("../data/flowering_time/floweringTimeMaizeDong2012.csv", stringsAsFactors = F, header = T, sep = ",") # small n=48 candidate set pulled from QTL and functional validation studies based on Arabidopsis flowering time pathway functions
flowering_pbe <- read.table("../data/flowering_time/TableS4_floweringTimeGlycerolipid.txt", stringsAsFactors = F, header = T, sep = "\t") # large n=903 set of candidate flowering time genes
flowering_gwas_male <- read.table("../data/flowering_time/male_flowering_time_Navarro2017.csv", stringsAsFactors = F, header = T, sep = ",")
flowering_gwas_female <- read.table("../data/flowering_time/female_flowering_time_Navarro2017.csv", stringsAsFactors = F, header = T, sep = ",")

# not too much overlap between sets, though more for the pathway sets
table(flowering_pathway$geneID %in% flowering_pbe$geneID)
table(flowering_pathway$geneID %in% flowering_pbe$geneID[flowering_pbe$sum >= 1])
table(flowering_gwas_male$ID %in% flowering_pbe$geneID)
table(flowering_gwas_female$ID %in% flowering_pbe$geneID)

# get conversion between coordinates downloaded from maizeGDB
g3_convert <- read.table("../data/refMaize/geneAnnotations/gene_model_xref_v3_from_maizeGDB.txt", 
                         stringsAsFactors = F, header = T, sep = "\t",
                         na.strings = c("NA", "<NA>", "")) %>%
  dplyr::select(c("v3_gene_model", "v3_start", "v3_end", "v3_chr", "v4_gene_model", "B73.Zm00001d.2."))
genes <- read.table("../data/refMaize/geneAnnotations/Zea_mays.B73_RefGen_v4.41.chr.genes.only.gff3", 
                    sep = "\t", 
                    stringsAsFactors = F)[ , c(1,4,5,9)]
colnames(genes) <- c("chr", "start", "end", "label")
genes$id <- sapply(substr(genes$label, 9, 100), function(x) strsplit(x, split = ";")[[1]][1])
genes$length = genes$end - genes$start + 1

# are my lists of genes in the .gff3? no. wrong gene name version.
table(flowering_gwas_female$ID %in% genes$id)
table(flowering_gwas_male$ID %in% genes$id)
table(flowering_pathway$geneID %in% genes$id)
table(flowering_pbe$geneID %in% genes$id) # none or on v4 gene sets
table(domestication$ID %in% genes$id) # none

# get conversion to other gene models
# oops! I'm not sure what's going on. Why doesn't my gff3 match v4? or match 
table(genes$id %in% g3_convert$B73.Zm00001d.2.)
table(genes$id %in% g3_convert$v4_gene_model) 
# same genes are in both gene models. but going the other way v4_gene_model has extra genes not in .gff3 
table(g3_convert$v4_gene_model == g3_convert$B73.Zm00001d.2.)
filter(g3_convert, v4_gene_model != B73.Zm00001d.2.) %>%
  mutate(match_found = v4_gene_model %in% genes$id) %>%
  dplyr::select(match_found) %>%
  table() # the ones that don't match B73.Zm00001d.2. aren't in the .gff3 file either
filter(g3_convert, v4_gene_model != B73.Zm00001d.2.) %>%
  dplyr::select(B73.Zm00001d.2.) %>%
  table() # no entry
# what is going on? gff3 has only half the genes in B73.Zm00001d.2.
table(unique(g3_convert$v4_gene_model) %in% genes$id)
table(unique(g3_convert$B73.Zm00001d.2.) %in% genes$id)

# convert to v4 list, id column must be ID
convert2v4 <- function(gene_set){
  id <- quo(id_col)
  left_join(gene_set, g3_convert, 
            by = c("ID"="v3_gene_model")) %>%
    filter(!is.na(B73.Zm00001d.2.)) %>% 
    dplyr::select(B73.Zm00001d.2.) %>%
    unique()
} 

flowering_pathway2 <- flowering_pathway %>%
  mutate(ID=geneID) %>%
  convert2v4()
flowering_pbe2 <- flowering_pbe %>%
  mutate(ID=geneID) %>%
  convert2v4()
flowering_gwas_male2 <- convert2v4(flowering_gwas_male)
flowering_gwas_female2 <- convert2v4(flowering_gwas_female)
domestication2 <- convert2v4(gene_set = domestication)


table(domestication_genes %in% g3_convert$B73.Zm00001d.2.)
table(domestication_genes %in% domestication2$B73.Zm00001d.2.)
table(unique(domestication2$B73.Zm00001d.2.) %in% domestication_genes) # I'm finding more genes. Shouldn't this be the same??
head(g3_convert[,c("v3_gene_model", "v4_gene_model")])
domestication_genes[duplicated(domestication_genes)]
domestication_genes[is.na(domestication_genes)]
domestication2[duplicated(domestication2$B73.Zm00001d.2),] # need to get rid of duplicates
# one example
domestication2[domestication2$B73.Zm00001d.2 == "Zm00001d015767", ]
genes[genes$id == "Zm00001d015767",]

# very few genes don't have hits v3 -> v4 gene models
genes$domestication <- genes$id %in% domestication_genes
genes$domestication2 <- genes$id %in% domestication2$B73.Zm00001d.2.
table(domestication_genes %in% genes$id) # about half the domestication genes can be found in the full genes list
table(domestication2$B73.Zm00001d.2. %in% genes$id) # slightly more included. more genes in general.
# but still only about half found in genes$id
# domestication2 has 89 more hits (though missing 10) .. going with that for now
table(genes[ , c("domestication", "domestication2")]) 
lapply(list(flowering_pathway2, flowering_pbe2, flowering_gwas_female2, flowering_gwas_male2),
       function(x) table(x$B73.Zm00001d.2. %in% genes$id)) # a bit more or less than half, depending on the set

# now I need to get an ancestry approximation for each gene
get_gene_anc <- function(anc_snps = d, gene, focal_col){
  anc_snps_chr <- dplyr::filter(anc_snps, chr == gene$chr)
  low <- ifelse(gene$start < anc_snps_chr$pos[1], 0, max(which(anc_snps_chr$pos < gene$start)))
  high <- ifelse(gene$end > anc_snps_chr$pos[nrow(anc_snps_chr)], nrow(anc_snps_chr), min(which(anc_snps_chr$pos > gene$end)))
  #return(data.frame(value = mean(anc_snps_chr[low:high, focal_col]), n = length(low:high)))
  return(mean(anc_snps_chr[low:high, focal_col]))
}
d %>%
  mutate(., start_bp = pos - 1) %>%
  mutate(., end_bp = pos) %>%
  dplyr::select(chr, start_bp, end_bp, meanAlpha_maize) %>%
  write.table(., "../data/domestication/meanAlpha_maize.bed", sep = "\t", col.names = F, row.names = F, quote = F)
d %>%
  mutate(., start_bp = pos - 1) %>%
  mutate(., end_bp = pos) %>%
  dplyr::select(chr, start_bp, end_bp, meanAlpha_mex) %>%
  write.table(., "../data/domestication/meanAlpha_mex.bed", sep = "\t", col.names = F, row.names = F, quote = F)
genes %>%
  filter(domestication) %>%
  dplyr::select(chr, start, end, id) %>%
  write.table(., "../data/domestication/domestication_genes.bed", sep = "\t", col.names = F, row.names = F, quote = F)
genes %>%
  filter(!domestication) %>%
  dplyr::select(chr, start, end, id) %>%
  write.table(., "../data/domestication/nondomestication_genes.bed", sep = "\t", col.names = F, row.names = F, quote = F)
genes %>%
  dplyr::select(chr, start, end, id) %>%
  write.table(., "../data/domestication/all_genes.bed", sep = "\t", col.names = F, row.names = F, quote = F)


# get mean ancestry for each gene: (VERY SLOW)
genes$meanAlpha_maize <- sapply(1:nrow(genes), function(i) get_gene_anc(anc_snps = d, gene = genes[i,], focal_col = "meanAlpha_maize"))
genes$meanAlpha_mex <- sapply(1:nrow(genes), function(i) get_gene_anc(anc_snps = d, gene = genes[i,], focal_col = "meanAlpha_mex"))

# after analysing with bedtools (see data/domestication/README.txt), reload into R
g <- read.table("../data/domestication/all_genes.meanAlpha_mex.autosomal.sorted.bed", header = F, stringsAsFactors = F, na.strings = ".")
colnames(g) <- c("chr", "start", "end", "gene", "meanAlpha_mex")
g10 <- read.table("../data/domestication/all_genes.meanAlpha_mex.10kb.autosomal.sorted.bed", header = F, stringsAsFactors = F, na.strings = ".")
colnames(g10) <- c("chr", "start", "end", "gene", "meanAlpha_mex")
g20 <- read.table("../data/domestication/all_genes.meanAlpha_mex.20kb.autosomal.sorted.bed", header = F, stringsAsFactors = F, na.strings = ".")
colnames(g20) <- c("chr", "start", "end", "gene", "meanAlpha_mex")

g_maize <- read.table("../data/domestication/all_genes.meanAlpha_maize.autosomal.sorted.bed", header = F, stringsAsFactors = F, na.strings = ".")
colnames(g_maize) <- c("chr", "start", "end", "gene", "meanAlpha_maize")
g10_maize <- read.table("../data/domestication/all_genes.meanAlpha_maize.10kb.autosomal.sorted.bed", header = F, stringsAsFactors = F, na.strings = ".")
colnames(g10_maize) <- c("chr", "start", "end", "gene", "meanAlpha_maize")
g20_maize <- read.table("../data/domestication/all_genes.meanAlpha_maize.20kb.autosomal.sorted.bed", header = F, stringsAsFactors = F, na.strings = ".")
colnames(g20_maize) <- c("chr", "start", "end", "gene", "meanAlpha_maize")

g10_maize %>%
  mutate(., domestication = gene %in% domestication_genes) %>%
  ggplot(aes(y = meanAlpha_maize, fill = domestication)) +
  geom_boxplot() +
  geom_abline(slope = 0, intercept = mean(d$meanAlpha_maize), 
              linetype = "dotted") +
  ggtitle("mean ancestry in maize pops near genes") +
  xlab("gene type") +
  ylab("mean mexicana ancestry in gene +/- 5kb")
ggsave("plots/mean_mex_anc_in_maize_near_domestication_genes_10kb.png",
       height = 8, width = 10, units = "in", device = "png")

g10 %>%
  mutate(., domestication = gene %in% domestication_genes) %>%
  ggplot(aes(y = meanAlpha_mex, fill = domestication)) +
  geom_boxplot() +
  geom_abline(slope = 0, intercept = mean(d$meanAlpha_mex), 
              linetype = "dotted") +
  ggtitle("mean ancestry in mexicana pops near genes") +
  xlab("gene type") +
  ylab("mean mexicana ancestry in gene +/- 5kb")
ggsave("plots/mean_mex_anc_in_mexicana_near_domestication_genes_10kb.png",
       height = 8, width = 10, units = "in", device = "png")


g10_maize %>%
  mutate(., domestication = gene %in% domestication_genes) %>%
  with(., summary(lm(meanAlpha_maize ~ domestication)))
g10 %>%
  mutate(., domestication = gene %in% domestication_genes) %>%
  with(., summary(lm(meanAlpha_mex ~ domestication)))
# domestication genes don't look any more likely than other genes
# to not have ancestry calls 
table(data.frame(domestication = g10$gene %in% domestication_genes, has_ancestry = is.na(g10$meanAlpha_mex)))

# flowering time
flowers <- c(flowering_gwas_female2$B73.Zm00001d.2.,
             flowering_gwas_male2$B73.Zm00001d.2.,
             flowering_pathway2$B73.Zm00001d.2.,
             flowering_pbe2$B73.Zm00001d.2.) %>%
  .[!duplicated(.)]
g %>%
  mutate(., domestication = gene %in% flowering_gwas_male2$B73.Zm00001d.2.) %>%
  ggplot(aes(y = meanAlpha_mex, fill = domestication)) +
  geom_boxplot() +
  geom_abline(slope = 0, intercept = mean(d$meanAlpha_maize), 
              linetype = "dotted") +
  ggtitle("mean ancestry in maize pops near genes") +
  xlab("gene type") +
  ylab("mean mexicana ancestry in gene +/- 5kb")



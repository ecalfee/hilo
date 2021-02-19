library(dplyr)
library(tidyr)
#setwd("~/Documents/gitErin/hilo")
# load variables from Snakefile
# unique flowering time candidate genes from Li et al 2016 supporting table S7
flowering_time_genes_Li_2016 = snakemake@input[["Li"]]
# flowering_time_genes_Li_2016 = "data/flowering_time/v2_unique_gene_names_Li_et_al_2016.csv"
# core 48 flowering time genes from known pathways
flowering_time_genes_Dong_2012 = snakemake@input[["Dong"]]
# flowering_time_genes_Dong_2012 = "data/flowering_time/floweringTimeMaizeDong2012.csv"
# all maize genes on v4 crossreferenced in other genome versions
cross_ref_file = snakemake@input[["cross_ref"]]
# cross_ref_file = "data/refMaize/geneAnnotations/gene_model_xref_v4_from_gramene.txt"
# bed file output
bed_file = snakemake@output[["bed"]]
# bed_file <- "ZAnc/results/flowering_time_genes_v4.bed"

# load data
# cross reference for all v4 genes
cross <- read.table(cross_ref_file, header = T,
                    stringsAsFactors = F, sep = "\t")
# core 48 flowering time genes 
core_flower <- read.table(flowering_time_genes_Dong_2012, sep = ",", 
                   header = T, stringsAsFactors = F) %>%
  rename(v2_gene_model = geneID) %>%
  mutate(v2_gene_model = trimws(v2_gene_model, "both")) %>% # trim any leading or trailing whitespace
  mutate(core_48 = T)
# larger list of 905 flowering time genes
larger_flower <- read.table(flowering_time_genes_Li_2016, sep = "\t", header = F,
                           stringsAsFactors = F) %>%
  rename(v2_gene_model = V1) %>% 
  mutate(v2_gene_model = trimws(v2_gene_model, "both")) # trim any leading or trailing whitespace

# all flowering time genes
flower <- larger_flower %>%
  # add all core flowering time genes (most overlap)
  full_join(., core_flower, by = "v2_gene_model") %>%
  mutate(core_48 = ifelse(is.na(core_48), F, core_48)) %>%
  # some candidate genes have no v2 gene name and therefore won't have v4 matches
  mutate(v2_gene_model = ifelse(v2_gene_model == "", NA, v2_gene_model)) %>% 
  # add v4 gene coordinates if found
  left_join(., cross, by = "v2_gene_model") %>%
  mutate(found_v4 = !is.na(v4_gene_model)) %>% # found gene on v4 genome
  mutate(on_chr_v4 = v4_chr %in% paste0("Chr", 1:10)) %>%
  mutate(ZmGene = ifelse(v2_gene_model == "GRMZM2G004483" & is.na(ZmGene), 
                         "ZmCCT9", ZmGene)) # identifying ZmCCT9 within the larger list
#table(duplicated(flower))
#dplyr::group_by(flower, found_v4, core_48) %>%
#  count()
# most, but not all candidate genes found in v4 gene model
# are localized to a chromosome (vs. unknown coordinates or on an unplaced contig)
#dplyr::filter(flower, found_v4) %>%
#  dplyr::group_by(., core_48, on_chr_v4) %>%
#    count()
#dplyr::filter(flower, found_v4 & !on_chr_v4) %>%
#  dplyr::arrange(., desc(core_48)) %>%
#  View(.) # look at exceptions

# note: 5 core flowering time genes from Dong 2012 aren't in the larger list 
# from Li 2016, or don't have a v2_gene_model name:
# the 3 with a name can be found on v4
# dplyr::filter(flower, !(v2_gene_model %in% larger_flower$v2_gene_model))

# make bed file with all flowering time candidate genes
# with coordinates found on v4 chromosomes 1-10
bed <- dplyr::filter(flower, on_chr_v4) %>%
  mutate(chr = substr(v4_chr, 4, 5), 
         start = v4_start, end = v4_end,
         gene_set = ifelse(core_48, "core_pathway", "all_candidate_genes")) %>%
  dplyr::select(chr, start, end, v4_gene_model, v2_gene_model, gene_set)
write.table(bed, file = bed_file, sep = "\t", quote = F, col.names = T, row.names = F)

#!/usr/bin/env Rscript
# Author: Erin Calfee 2018
# This script assigns the labels 'maize' 'mexicana' 'parviglumis'
# to the ancestries that come out of NGSadmix, based
# on which group has the most of that ancestry
library(dplyr)
library(tidyr)

# take in command line arguments
args = commandArgs(trailingOnly=TRUE)

# to run: Rscript define_K_ancestries.R pass2_alloMAIZE_PalmarChico 1 1cM 1 3

#PREFIX <- "pass2_alloMAIZE_PalmarChico"
PREFIX = args[1]
#r = 1 # recombination bin
r = as.numeric(args[2])
#WIND = "1cM"
WIND = args[3]
#BOOT = 0 # bootstrap
BOOT = args[4]
#K = 3
K = as.numeric(args[5])
PREFIX_METRICS <- "hilo_alloMAIZE_MAIZE4LOW"

ancestries <- c("maize", "mexicana", "parviglumis")[1:K]
file_prefix <- paste0("results/bootstrap/windows_", WIND,
                      "/r5_recomb", r, "/", PREFIX, "/K", K)
admix <- read.table(paste0(file_prefix, "/boot", BOOT, ".qopt"))
colnames(admix) <- paste0("anc", 1:K) #c("anc1", "anc2")
file_out <- paste0(file_prefix, "/boot", BOOT, ".anc")

# get labels
IDs <- data.frame(ID = read.table(paste0("../samples/", PREFIX, "_IDs.list"), header = F, stringsAsFactors = F)$V1, stringsAsFactors = F)
metrics <- read.table(paste0("../filtered_bams/metrics/", PREFIX_METRICS, ".flagstat.total"), header = F, stringsAsFactors = F)
colnames(metrics) <- c("ID", "total_reads_pass")
hilo <- read.table("../samples/hilo_meta.txt", stringsAsFactors = F, header = T, sep = "\t")
landraces <- read.table("../samples/alloMAIZE_meta.txt", stringsAsFactors = F, header = T, sep = "\t")
parviglumis <- read.table("../samples/PalmarChico_meta.txt", stringsAsFactors = F, header = T, sep = "\t")

# a rough coverage estimate is number of reads passing filtering * 150bp/read
# divided by total area they could map to in maize reference genome v4

# size of reference genome reads are mapped to
ref_genome_size <- sum(read.table("../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.fai", 
                                  stringsAsFactors = F)$V2) # V2 is the size of each chromosome mapped to,
# including parts I won't analyze on the Pt and Mt and scaffolds not assigned to chromosomes (but reads would pass Q filters there too)
metrics$est_coverage = round(metrics$total_reads_pass*150/ref_genome_size, 4)

# combine sample meta data
meta <- bind_rows(hilo, landraces, parviglumis) %>%
  left_join(., metrics, by = "ID") %>%
  mutate(., group = paste(symp_allo, zea, sep = "_"))

# join bams and admix by position (CAUTION - bam list order and admix results MUST MATCH!)
d <- bind_cols(IDs, admix)  %>%
  left_join(., meta, by = "ID") %>%
  arrange(., popN) %>%
  arrange(., zea) %>%
  arrange(., symp_allo)

min_coverage = 0.1

which_anc <- data.frame(ancestry = colnames(admix),
                   ancestry_label = 
                     sapply(colnames(admix), function(x) 
                       names(which.max( # exlcude lowest coverage samples
                         tapply(d[d$zea == "parviglumis" | d$est_coverage > min_coverage, x], 
                                d$zea[d$zea == "parviglumis" | d$est_coverage > min_coverage], 
                                mean)))),
                   stringsAsFactors = F)

d %>% 
  tidyr::gather(., "ancestry", "p", colnames(admix)) %>%
  left_join(., which_anc, by = "ancestry") %>%
  dplyr::select(-ancestry) %>%
  tidyr::spread(., ancestry_label, p) %>%
  dplyr::select(ID, ancestries) %>%
  write.table(., file_out,
              quote = F,
              sep = "\t", col.names = T, row.names = F)
  



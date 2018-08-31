#!/usr/bin/env Rscript

# output:
# this script makes a file with allopatric allele counts
# and SNP information, and genetic map distance between markers
# this output is used by make_input_ancestry_hmm.R, which adds population-specific
# read counts for the admixed population being run

# inputs:
# random seed for sampling 
# for high-coverage allopatric maize it uses allele counts
# from files .frq.counts files generated in plink
# for sympatric maize it samples 1 read each per individual
# from a read count file e.g. hilo_N_chr1.csv produced by ASEReadCounter (GATK)
# SNP information is obtained from .var.sites file
# and map information (distance between SNPs in Morgans) is obtained from .distM

library(dplyr)
# to run:
# Rscript make_allo_counts_ancestry_hmm.R 4 ../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM

print(getwd()) # print current directory
# set random seed
rseed = 13248934
set.seed(rseed)

# arguments
args = commandArgs(trailingOnly=TRUE)
# chromosome #
#i = 4
i = as.integer(args[1])
# path to input data
#path = "../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM"
path = args[2]

# helper function
# takes in GATK file and returns counts data frame with 'chr' 'pos' 
# 'refCount' and 'altCount' columns after sampling 1 read per individual
# with either 0 0, 1 0 or 0 1 as values for refCount and altCount at each SNP
sample1_read = function(GATK_counts_file){
  counts = read.table(GATK_counts_file, header = T, stringsAsFactors = F, sep = "\t") %>%
      rename(., chr = contig) %>%
      rename(., ref = refAllele) %>%
      rename(., alt = altAllele) %>%
      select(., chr, position, ref, alt, refCount, altCount)
  tot = counts$refCount + counts$altCount
  # if any reads are observed (tot > 0), sample 1 read. Otherwise count zero reads.
  new_altCount = sapply(1:length(tot), function(x) ifelse(tot[x] == 0, 0, rbinom(n = 1,
                                            size = 1, # always sample 1 read
                                            prob = counts$altCount[x]/tot[x])))
  new_refCount = ifelse(tot == 0, 0, 1 - new_altCount)
  new_counts = counts %>%
    mutate(., refCount = new_refCount) %>%
    mutate(., altCount = new_altCount)
  if (sum(is.na(new_counts)) > 0){ # there should not be any NA values at this point
    print(new_counts[!complete.cases(new_counts), ])
    stop("there are some NAs in counts")
  }
  return(new_counts)
}

# input data
SNPs = read.table(paste0(path, "/chr", i, ".var.sites"),
                  header = F, stringsAsFactors = F, sep = "\t")
colnames(SNPs) = c("chr", "position", "ref", "alt")
map_pos = read.table(paste0(path, "/chr", i, ".distM"), 
                     header = F, stringsAsFactors = F)
colnames(map_pos) = c("distM")

# maize counts
# maize vcf dropped some positions (presumably where no genotypes could be called)
# therefore maize_pos is shorter than SNPs..which I deal with below and recalculate
# genetic distances between positions kept
maize_pos <- read.table(paste0(path, "/maize.allo.4Low16_chr", i, ".vcf"),
                        header = F, sep = "", stringsAsFactors = F)[, c(1,2,4,5)]
colnames(maize_pos) = c("chr", "position", "ref", "alt")
maize_counts <- read.table(paste0(path, "/maize.allo.4Low16_chr", i, ".frq.count"),
                           header = T, sep = "", stringsAsFactors = F) %>%
  mutate(., refCount = C2) %>%
  mutate(., altCount = C1) %>%
  select(., refCount, altCount) %>%
  bind_cols(maize_pos, .)

# mexicana counts
# find allopatric individuals to include
pass1 <- read.table("../data/pass1_ids.txt", stringsAsFactors = F, 
                             header = T, sep = "\t")
# allopatric mexicana includes Xochimilco (pop # 35) for now
mex_ids <- filter(pass1, zea == "mexicana" & 
                    (symp_allo == "allopatric" | LOCALITY == "Xochimilco"))$n
altCountsAll = rep(0, nrow(SNPs))
refCountsAll = rep(0, nrow(SNPs))
for (id in mex_ids){
  newCounts = left_join(SNPs, sample1_read(
      GATK_counts_file = paste0(path, "/hilo_", id, "_chr", i, ".csv")),
    by = c("chr", "position", "ref", "alt"))
  newCounts[is.na(newCounts)] <- 0 # NA values become 0 ref and 0 alt counts (NA b/c SNP is missing for individual -- no coverage)
  altCountsAll = altCountsAll + newCounts$altCount
  refCountsAll =  refCountsAll + newCounts$refCount 
}  
mex_counts = data.frame(altCount = altCountsAll, refCount = refCountsAll,
                        stringsAsFactors = F)
d1 = SNPs %>%
  left_join(., maize_counts, by = c("chr", "position", "ref", "alt")) %>%
  bind_cols(., mex_counts[ , c("refCount", "altCount")], map_pos) %>%
  mutate(., posM = cumsum(distM)) %>%
  filter(., !is.na(refCount)) %>%
  filter(., !(refCount == 0 & altCount == 0)) %>% # no genotypes called in maize 
  filter(., !(refCount1 == 0 & altCount1 == 0)) %>% # no informative sites from mexicana -- need to check upstream why these aren't being pre-filtered out
  mutate(., distM_new = c(1, diff(posM))) %>% # after filtering out uninformative sites, get new genetic distances between remaining sites
  select(., c("chr", "position", "ref", "alt", "refCount", "altCount", "refCount1", "altCount1", "distM_new")) %>%
  rename(., refCount_maize = refCount) %>%
  rename(., altCount_maize = altCount) %>%
  rename(., refCount_mex = refCount1) %>%
  rename(., altCount_mex = altCount1)
#head(with(d1, which(distM_new != distM))) most don't match due to small rounding errors
# this is ok for preliminary analysis now, but I need a better work-around for future
# picking SNPs where they won't drop out when genotyping maize AND/OR
# using plink to upload map distances at the end after filtering
# for now the error is only 15 or 16 digits in .. so I will format to include 10 digits only

# write output to be used by make_input_ancestry_hmm.R
write.table(format(d1, digits = 10), paste0(path, "/allo_counts_chr", i, ".txt"),
            row.names = F, col.names = T, quote = F)
warnings()

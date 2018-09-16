#!/usr/bin/env Rscript

# output:
# this script makes a file with allopatric allele counts
# and SNP information, and genetic map distance between markers
# this output is used by make_input_ancestry_hmm.R, which adds population-specific
# read counts for the admixed population being run

# inputs:
# (1) random seed for sampling
# (2) for high-coverage allopatric maize it uses allele counts
# from files .frq.counts files generated in plink
# (3) for allopatric reference maize it samples 1 read each per individual
# from a major/minor allele read count file produced by running
# countReadsACGT.sh then countReadsMajorMinor.sh for each reference individual
# (4) SNP information is obtained from .var.sites file
# (5) and map information (distance between SNPs in Morgans) is obtained from .distM

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
# takes in individual counts file and samples 1 read per individual
# with either 0 0, 1 0 or 0 1 as values for n_major and n_minor at each SNP
sample1_read = function(ind_counts_file){
  counts = read.table(ind_counts_file, header = F, stringsAsFactors = F, sep = "\t")
  colnames(counts) = c(n_major, n_minor)
  tot = counts$n_major + counts$n_minor
  # if any reads are observed (tot > 0), sample 1 read. Otherwise count zero reads.
  new_n_minor = sapply(1:length(tot), function(x) ifelse(tot[x] == 0, 0, rbinom(n = 1,
                                            size = 1, # always sample 1 read
                                            prob = counts$n_minor[x]/tot[x])))
  new_n_major = ifelse(tot == 0, 0, 1 - new_n_minor)
  new_counts = data.frame(n_major = new_n_major, n_minor = new_n_minor)
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
colnames(maize_pos) = c("chr", "position", "major", "minor")
maize_counts <- read.table(paste0(path, "/maize.allo.4Low16_chr", i, ".frq.count"),
                           header = T, sep = "", stringsAsFactors = F) %>%
  mutate(., n_major_maize = C2) %>%
  mutate(., n_minor_maize = C1) %>%
  select(., n_major_maize, n_minor_maize) %>%
  bind_cols(maize_pos, .)

# mexicana counts
# find allopatric individuals to include
pass1 <- read.table("../data/pass1_ids.txt", stringsAsFactors = F,
                             header = T, sep = "\t")
# allopatric mexicana includes Xochimilco (pop # 35) for now
mex_ids <- filter(pass1, zea == "mexicana" &
                    (symp_allo == "allopatric" | LOCALITY == "Xochimilco"))$n
mex_counts = data.frame(n_major_mex = rep(0, nrow(SNPs)), n_minor_mex = rep(0, nrow(SNPs)))
for (id in mex_ids){
  newCounts = sample1_read(
      ind_counts_file = paste0(path, "/countsMajMin/chr", i, "/hilo_", id, "counts.txt"))
  # update count totals
  mex_counts$n_minor_mex = mex_counts$n_minor_mex + newCounts$n_minor
  mex_counts$n_major_mex = mex_counts$n_major_mex + newCounts$n_major
}

d = SNPs %>%
   left_join(., maize_counts, by = c("chr", "position", "major", "minor") %>%
   bind_cols(., mex_counts, map_pos) %>%
   #mutate(., posM = cumsum(distM)) %>% # code to exclude sites with no genotypes called in maize
   #filter(., !is.na(n_major_maize)) %>% # no genotypes called in maize
   #mutate(., distM_new = c(1, diff(posM))) %>% # after filtering out uninformative sites (no maize genotypes), get new genetic distances between remaining sites
   #mutate(., distM = distM_new) %>%
   select(., c("chr", "position", "major", "minor", "n_major_maize", "n_minor_maize",
   "n_major_mex", "n_minor_mex", "distM"))

print(paste0("maize has no genotype called for ",
sum(is.na(d$n_major_maize), "SNPs")) # print SNPs that have no genotypes..despite min # ind's with data
print(head(d[is.na(d$n_major_maize, ])) # these SNPs are not excluded
# write output to be used by make_input_ancestry_hmm.R (with no headers)
write.table(format(d, digits = 10), paste0(path, "/allo_counts_chr", i, ".txt"),
    row.names = F, col.names = T, quote = F) # prints column names even though end output won't
warnings()

#!/usr/bin/env Rscript
source("calcMapLinearApprox.R") # linear interpolation functions

# this script makes windows of a specific cM width
# around each SNP position given

print("current working directory:")
print(getwd()) # print current directory

# to run:
# Rscript makeWindowsCM.R CM_WINDOW SNP_FILE_IN BED_FILE_OUT

# input: 
args = commandArgs(trailingOnly=TRUE)
# (1) a cM distance for desired width of each window
CM_WINDOW = as.numeric(args[1])
#CM_WINDOW = 0.1
# (2) a tab-deliminated SNP position file (no header) with chromosome and bp position as 1st and 2nd columns
SNP_FILE_IN = args[2]
#SNP_FILE_IN = "../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM/chr6.var.sites"
BED_FILE_OUT= args[3]
#BED_FILE_OUT = paste0(tools::file_path_sans_ext(SNP_FILE_IN), "_", CM_WINDOW, "cM_windows.bed")

# output: a bed file BED_FILE_OUT
# with each bp interval and the total cM width of the window
# (end cases are taken to the end of the chromosome)
# cM distances are calculated using the Ogut 2015 maize genetic map
# extended to the ends of each chromosome
# note: I follow convention making bed start positions
# zero-based and end positions 1-based
# as in positions 1-5 would be indicated 0-5 in output.bed
# GFF3 files have all 1-based positions (but bedtools deals with this)

# load chromosome start and end positions from extended map
chr_starts = rmapEXT[rmapEXT$pos_bp == 1, ]
chr_ends = rmapEXT[sapply(unique(rmapEXT$chr),
                          function(i) which(rmapEXT$chr == i & 
                                              rmapEXT$pos_bp == max(rmapEXT[rmapEXT$chr == i, "pos_bp"]))), ]

SNPs = read.table(SNP_FILE_IN,
           stringsAsFactors = F, header = F, sep = "\t")[ , 1:2]
colnames(SNPs) = c("chrom", "pos_bp")

# get center, left and right endpoints for windows
# in cM
center_cm = apply(SNPs, 1, function(row)
  bp2cM(chrom = row["chrom"], pos_bp = row["pos_bp"]))
left_cm = center_cm - CM_WINDOW/2
right_cm = center_cm + CM_WINDOW/2

# translate to cM to bp
left_bp = mapply(FUN = function(i, p) 
  cM2bp(chrom = i, pos_cM = p),
                 SNPs$chrom, left_cm)
right_bp = mapply(FUN = function(i, p) 
  cM2bp(chrom = i, pos_cM = p),
  SNPs$chrom, right_cm)

# edge effects -- fix NA's (outside of chromosome boundaries)
# low outliers take first position on chromosome
low = which(is.na(left_bp))
for (i in low){
  left_bp[i] <- 1
  left_cm[i] <- chr_starts[chr_starts$chr == SNPs$chrom[i], "pos_cM"]
}
# high outliers take last position on chromosome
high = which(is.na(right_bp))
for (i in high){
  right_bp[i] <- chr_ends[chr_ends$chr == SNPs$chrom[i], "pos_bp"]
  right_cm[i] <- chr_ends[chr_ends$chr == SNPs$chrom[i], "pos_cM"]
}
# round to 6 decimal places -- realistically we don't have lower precision
width_cm = round(right_cm - left_cm, 6)

# write a bed file

write.table(data.frame(chrom = SNPs$chrom,
                       chromStart = left_bp - 1, # standard bed files start is index 0 
                       chromEnd = right_bp, # and end is index 1
                       name = width_cm), # name field saves cM width of window (could be short if near edge of chromosome)
            file = BED_FILE_OUT,
            sep = "\t", quote = F, 
            col.names = T, row.names = F)


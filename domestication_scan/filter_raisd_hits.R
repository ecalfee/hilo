#!/usr/bin/env Rscript
# working directory is hilo/
library(GenomicRanges)
library(data.table)
library(dplyr)
library(seqinr)
# script filters raisd hits for too many Ns (missing data in reference genome)

# input files
bed_in <- snakemake@input[["bed"]]
# bed_in <- "domestication_scan/results/raisdHits.bed"
fa_in <- snakemake@input[["fa"]]
# fa_in <- "domestication_scan/results/raisdHits.fa"

# parameters
maxN <- as.numeric(snakemake@params[["maxN"]])
# maxN <- 0.5 # proportion to use as cutoff for maximum Ns

# output files
bed_keep <- snakemake@output[["bed_keep"]]
# bed_keep <- "domestication_scan/results/raisdHits_keep.bed"
bed_excl <- snakemake@output[["bed_excl"]]
# bed_excl <- "domestication_scan/results/raisdHits_exclude.bed"


# check to see which of the outliers are driven by N's

# read in the bed file:
hits <- data.table::fread(bed_in, data.table = FALSE)

# read in the sequences of the regions:
fas <- read.fasta(fa_in)

# calculate percent of bp with missing data in the reference
Ns <- sapply(fas, function(x) sum(x == 'n'))
lengths <- sapply(fas, function(x) length(x))
percents <- Ns/lengths

# filter out raisd hits exceeding maximum missing data (likely false-positives)
hits_keep <- hits[which(percents < maxN), ] 
hits_excl <- hits[which(percents >= maxN), ]

# sum(hits_keep$V3 - hits_keep$V2)/10^9/2.1 # covers about 2% of the ref genome
print("percentage of 2.1Gb genome:")
print(sum(hits_keep$V3 - hits_keep$V2)/10^9/2.1)

write.table(hits_keep, file = bed_keep,
            quote = FALSE, sep = '\t',
            row.names = FALSE, col.names = FALSE)

write.table(hits_excl, file = bed_excl,
            quote = FALSE, sep = '\t',
            row.names = FALSE, col.names = FALSE)

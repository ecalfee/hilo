#!/usr/bin/env Rscript
# working directory is hilo/
library(GenomicRanges)
library(data.table)
library(dplyr)
# script creates bed file of top X% raisd hits (as contiguous regions) from raw raisd output

# input files
raisd_input <- snakemake@input[["raisd"]]
# raisd_input <- "domestication_scan/results/RAiSD_Output2.Raisd"
# parameters
top_x <- as.numeric(snakemake@params[["top_x"]])
# top_x <- 0.0005 # proportion to use as cutoff
# output files
bed_out <- snakemake@output[["bed"]]
# bed_out <- "domestication_scan/results/raisdHits.bed"


# load in raw Raisd results
raisd <- data.table::fread(raisd_input, data.table = FALSE)
colnames(raisd) <- c('seqnames','testsnp','start','end','a','b','c','score')
raisd$seqnames <- paste('chr', raisd$seqnames, sep="")

# Sort
raisdSort <- raisd[order(raisd$score, decreasing=TRUE), ]
# Rough attempt to get top X%
raisdHits <- raisdSort[1:round(nrow(raisd)*top_x), ] #testing because the final length will depend on how much overlap is here

# reduce hits with overlapping windows to contiguous ranges
grRaisd <- makeGRangesFromDataFrame(raisdHits, keep.extra.columns=TRUE) %>%
  reduce()

# write a bed file of grRaisd (will later remove any hits that are mostly Ns)
rdf <- data.frame(grRaisd)
rdf$seqnames <- sapply(rdf$seqnames, substr, 4, 6) # chr2 -> 2

# sum(rdf$width/10^9)/2.1 # top_x = 0.0005 covers about 3.6% of the ref genome before filtering out high N (missing data) windows
print("proportion of 2.1Gb genome:")
print(sum(rdf$width/10^9)/2.1)

write.table(rdf[ , c("seqnames", "start", "end")], file = bed_out,
            quote = FALSE, sep = '\t',
            row.names = FALSE, col.names = FALSE)

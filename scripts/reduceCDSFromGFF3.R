# installing some new packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("BiocInstaller", version = "3.8")
#library(BiocManager)
#BiocManager::install("IRanges")
library(IRanges)
library(dplyr)
d <- read.csv("../data/refMaize/geneAnnotations/Zea_mays.B73_RefGen_v4.41.chr.gff3",
                header = F, stringsAsFactors = F, skip = 20,
                sep = "\t")
colnames(d) <- c("chrom", "gramene", "type", "start_bp", "end_bp", "V6", "V7", "V8", "ID")

cds <- d %>%
  select(chrom, type, start_bp, end_bp) %>%
  filter(type == "CDS") # only CoDing Sequences, not whole exons

# many are still overlapping or even start at the same bp
# I need to find a way to either expand over each bp they cover and remove duplicates
# (maybe one chrom at a time to save space?)
# or to otherwise group by overlapping sequences and pick 
# the min start and max end between them
# looks like IRanges() could help me
#?IRanges()
# IRanges objects store start end and width or regions
# and reduce() collapses overlapping ranges into the largest contiguous ranges
# below I create a list of reduced irange objects, one per chromosome
chromosomes = c(1:10, "Mt", "Pt")
ir = lapply(chromosomes, # reduce each chrom separately
            function(chrom)
              IRanges::reduce(IRanges(start = cds[cds$chrom == chrom, "start_bp"], 
                    end = cds[cds$chrom == chrom, "end_bp"])))
# how many CDS per chromosome?
lapply(ir, length)
# print a file of coding ranges for each chromosome
lapply(1:12, 
       function(i) write.table(ir[[i]], 
                   paste0("../data/refMaize/geneAnnotations/CDS_", chromosomes[i] , ".txt"),
                   quote = F,
                   col.names = T,
                   row.names = F)
       )


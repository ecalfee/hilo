#!/usr/bin/env Rscript

# this script converts var sites used in hmm ancestry inference
# into tracts around those sites

# to run:
# Rscript ancestry_posterior_to_tracts.R ../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM

# arguments
args = commandArgs(trailingOnly=TRUE)
# directory with chri.var.sites files
input_dir = args[1]

# helper function makes start and end positions of tracts around a SNP site
# as the bp position halfway between that marker and the next closest marker.
# or, if marker has no nearest neighbor (at either end of chromosome),
# tract ends at the marker on one side
get_start_end <- function(chrom, dir){
  sites <- read.table(paste0(dir, "/chr", chrom, ".var.sites"),
                      stringsAsFactors = F, header = F)
  colnames(sites) <- c("chrom", "pos", "allele1", "allele2")
  sites$start = floor(sites$pos - diff(c(sites$pos[1]-1, sites$pos), lag = 1)/2)
  sites$end = floor(sites$pos + diff(c(sites$pos, sites$pos[length(sites$pos)]), lag = 1)/2)
  return(sites[, c("chrom", "start", "end")])
}
# find start and end positions for all sites across 10 chromosomes
all_sites <- do.call(rbind,
                     lapply(1:10, function(i)
                       get_start_end(chrom = i, dir = input_dir)))
# create output directory
dir.create(file.path("results", "input"))
# write sites to output file
write.table(all_sites, file.path("results", "input", "var.sites.bed"),
            col.names = F, row.names = F, quote = F)


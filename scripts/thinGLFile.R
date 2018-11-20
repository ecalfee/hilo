#!/usr/bin/env Rscript

# this script takes in a GL file in DIR_GL_IN
# with all samples and all SNPs for a region of the genome,
# and outputs a 'thinned' GL file to DIR_OUT
# with the same samples but only a subset of the SNPs specified by
# chromosome-level var.sites file in DIR_THINNED_SITES

# to run:
# Rscript thinGLFile.R REGION DIR_GL_IN DIR_THINNED_SITES DIR_OUT

# helper function
# split marker chr_positionBP into chr and pos_bp; return chromosome
chr_from_marker = function(marker){
  chrom = as.numeric(strsplit(marker, "_")[[1]])[1]
  return(chrom)
}


print(getwd()) # print current directory

# arguments
args = commandArgs(trailingOnly=TRUE)


# set variables
REGION = as.integer(args[1])
#REGION = 406
DIR_GL_IN = args[2]
#DIR_GL_IN = "../data/geno_lik/merged_pass1_all_alloMaize4Low_16/allVar_depthFilt"
FILE_GL_IN = paste0(DIR_GL_IN, "/region_", REGION, ".beagle.gz")
DIR_THINNED_SITES = args[3]
#DIR_THINNED_SITES = "../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedPCA"
DIR_OUT = args[4]
FILE_OUT = paste0(DIR_OUT, "/region_", REGION, ".beagle.gz")

# get genotype likelihoods for current region
region_gl <- read.table(gzfile(FILE_GL_IN),
                        header = T, stringsAsFactors = F)

# which chromosomes are in my gl file?
gl_chroms = unique(sapply(region_gl$marker,
                          function(i) chr_from_marker(i)))

# get list of thinned sites for all present chromosomes
thin_sites <- do.call(rbind,
                      lapply(gl_chroms, function(i)
                        read.table(paste0(DIR_THINNED_SITES, "/chr", i, ".var.sites"),
                         stringsAsFactors = F, header = F)))
colnames(thin_sites) <- c("chr", "pos_bp", "Major", "Minor")

# merge chr and pos_bp to create a marker column
# that matches the gl file marker column
thin_sites$marker = with(thin_sites,
                         paste(chr, pos_bp, sep = "_"))

# keep only gl's for SNPs in 'thin_sites'
thin_gl = region_gl[region_gl$marker %in% thin_sites$marker,]

# write thinned gl to gzipped output file
write.table(thin_gl, gzfile(FILE_OUT),
            quote = F, col.names = T, row.names = F)

# print any warnings
warnings()

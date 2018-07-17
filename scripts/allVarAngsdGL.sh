#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J allGLAngsd
#SBATCH -o /home/ecalfee/hilo/slurm-log/allGLAngsd_%j_%A_%a.out
#SBATCH -t 36:00:00
#SBATCH --mem=40G
#SBATCH -n 4
#SBATCH --array=0-46
#SBATCH --export=DIR_OUT=geno_lik/pass1/allVar,BAM_IN=pass1_bam.all.list,CHECK_BAM_HEADERS=1

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load angsd
module load angsd

# make directory to store output (if doesn't yet exist)
mkdir -p $DIR_OUT

# apply filtering with SAMtools & PICARD
echo "calling variants and GL using ANGSD on BAMS for hilo genomic regions "$SLURM_ARRAY_TASK_ID
# steps:
# (0) Start with filtered BAM files and reference genome
# (1) For each chromosomal region individually, find variant sites
angsd -out $DIR_OUT/region_$SLURM_ARRAY_TASK_ID \
-r $(cat refMaize/divide_50Mb/region_$SLURM_ARRAY_TASK_ID.txt) \
-rf refMaize/divide_50Mb/region_$SLURM_ARRAY_TASK_ID.txt \
-ref /group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa \
-checkBamHeaders $CHECK_BAM_HEADERS \
-bam $BAM_IN \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doMajorMinor 2 \
-doCounts 1 -minMaf 0.05 -doMaf 8 \
-GL 1 -doGlf 2 \
-minInd 40 \
-P 4


# settings:
# -checkBamHeader 1 means it checks all chrom/etc. match up between bam files included; 0 won't do the check. 
# I will set to zero for some cases, e.g. b/c hilo bams have mt chromosome but allopatric maize don't; I manually checked other chroms are compatible
# -r specifies which region to work on; -rf is the regions file
# -remove_bads removes reads with flags like duplicates
# -doMajorMinor 2: infer major and minor from allele counts
# -bam list of bams to include (all newly sequenced allopatric mex. and sympatric mexicana & maize pops)
# -GL 1: use samtools genotype likelihood method
# -doGlf 2: output beagle likelihood file
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# we use AlleleCounts method for MAF (-doMaf 8) with -doCounts
# which doesn't consider base quality and ignores all reads with non-target alleles
# but also doesn't rely on HW equlibrium assumptions to get ML allele freq from genotype freqs
# -minMaf x: and then do a cutoff to only include variant sites with >x minor allele freq.
# -minInd N: only keep sites with information (at least one read) from N individuals
# -P 4 means use 4 threads/nodes for each angsd task (here task=chromosome; then merges threads within-chrom)
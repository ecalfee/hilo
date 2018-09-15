#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J varAngsdReg
#SBATCH -o /home/ecalfee/hilo/slurm-log/varAngsdRegions_%A_%a.out
#SBATCH -t 10:00:00
#SBATCH --mem=12G
#SBATCH -n 4
#SBATCH --export=OUT_DIR=var_sites/pass1_alloMaize4Low,BAM_LIST=merged_bam.pass1_all.alloMaize4Low_16.list
#SBATCH --array=0-46

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load angsd -- don't load; use local copy v9.20
# module load angsd

# make directory to store output (if doesn't yet exist)
mkdir -p $OUT_DIR

# get variant sites and their MAF
echo "calling variants using ANGSD on BAMS for hilo region "$SLURM_ARRAY_TASK_ID
# steps:
# (0) Start with filtered BAM files and reference genome
# (1) For each chromosome individually, find variant sites
angsd -out $OUT_DIR/region_$SLURM_ARRAY_TASK_ID \
-rf refMaize/divide_50Mb/region_$SLURM_ARRAY_TASK_ID.txt \
-doMajorMinor 2 \
-bam $BAM_LIST \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doCounts 1 -minMaf 0.05 -doMaf 8 \
-minInd 50 \
-P 4


# settings:
# -r specifies which region to work on
# -remove_bads removes reads with flags like duplicates
# -doMajorMinor 4: pre-specify major allele from reference genome and infer minor allele from genotype likelihood
# -bam list of bams to include (all newly sequenced allopatric mex. and sympatric mexicana & maize pops)
# -GL 1: use samtools genotype likelihood method
# -doGlf 2: output beagle likelihood file
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# -doMaf 2 : output minor allele freq
# an alternative would be to use AlleleCounts method (-doMaf 8) with -doCounts
# which doesn't consider base quality and ignores all reads with non-target alleles
# but also doesn't rely on HW equlibrium assumptions to get ML allele freq from genotype freqs
# -minMaf x: and then do a cutoff to only include variant sites with >x diff. from reference allele
# -minInd N: only keep sites with information (at least one read) from N individuals
# -P 4 means use 4 threads/nodes for each angsd task (here task=chromosome; then merges threads within-chrom)

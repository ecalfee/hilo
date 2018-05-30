#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J varAngsd
#SBATCH -o /home/ecalfee/hilo/slurm-log/varAngsd_%j.out
#SBATCH -t 24:00:00
#SBATCH -x bigmem1
#SBATCH -mem 30G
#SBATCH -n 10

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load angsd
module load angsd

# make directory to store output (if doesn't yet exist)
mkdir -p var_sites

# apply filtering with SAMtools & PICARD
echo "calling variants using ANGSD on BAMS for hilo $SLURM_ARRAY_TASK_ID"
# steps:
# (0) Start with filtered BAM files and reference genome
# (1) For each chromosome individually, find variant sites

angsd -out var_sites/pass1 \
-r 1:10000-20000 \
-doMajorMinor 4 -ref /group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa \
-bam data/pass1_bam.all.list \
-GL 1 -doGlf 2 \
-baq 1 -minMapQ 30 -minQ 20 \
-minMaf 0.1 -doMaf 2 \
-minInd 50 \
-P 10

# settings:
# -r specifies which region to work on
# -doMajorMinor 4: pre-specify major allele from reference genome and infer minor allele from genotype likelihood
# -bam list of bams to include (all newly sequenced allopatric mex. and sympatric mexicana & maize pops)
# -GL 1: use samtools genotype likelihood method
# -doGlf 2: output beagle likelihood file
# -baq 1 -minMapQ 30 -minQ 20: calculate BAQ adjusted quality score and filter out sites with low mapping quality or base/BAQ quality
# -doMaf 2 : output minor allele freq
# an alternative would be to use AlleleCounts method (-doMaf 8)
# -minMaf 0.1: and then do a cutoff to only include variant sites with >10% diff. from reference allele
# -minInd 50: only keep sites with information (at least one read) from 50 individuals
# -P 10 sets 10 threads for use by ANGSD

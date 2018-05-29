#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J varAngsd
#SBATCH -o /home/ecalfee/hilo/slurm-log/varAngsd_%j_%A_%a.out
#SBATCH -t 24:00:00
#SBATCH -x bigmem1
#SBATCH --array=1-10
#SBATCH -mem 30G

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load angsd
module load angsd

# make directory to store output (if doesn't yet exist)
mkdir -p var_sites
# make a local ‘scratch’ directory for temporary files (@ end check that it’s empty)
mkdir -p /scratch/ecalfee/hilo_chr$SLURM_ARRAY_TASK_ID  # temporary sort files will be written to local node

# apply filtering with SAMtools & PICARD
echo "calling variants using ANGSD on BAMS for hilo $SLURM_ARRAY_TASK_ID"
# steps:
# (0) Start with filtered BAM files and reference genome
# (1) For each chromosome individually, find variant sites

angsd -out var_sites/chr$SLURM_ARRAY_TASK_ID \
-r $SLURM_ARRAY_TASK_ID \
-doMajorMinor 4 -ref /group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa \
-bam data/pass1_bam.all.list \
-GL 1 -doGlf 2 \
-baq 1 -minMapQ 30 -minQ 20 \
-doMaf 4 -minMaf 0.1 \
-rmTriallelic 1e-5 \
-minInd 50

# settings:
# -r specifies which region to work on
# -doMajorMinor 4: pre-specify major allele from reference genome and infer minor allele from genotype likelihood
# -bam list of bams to include (all newly sequenced allopatric mex. and sympatric mexicana & maize pops)
# -GL 1: use samtools genotype likelihood method
# -doGlf 2: output beagle likelihood file
# -baq 1 -minMapQ 30 -minQ 20: calculate BAQ adjusted quality score and filter out sites with low mapping quality or base/BAQ quality
# -doMaf 4 : use frequency from genotype probabilities to calculate minor allele freq
# an alternative would be to use AlleleCounts method (-doMaf 8)
# -minMaf 0.1: and then do a cutoff to only include variant sites with >10% diff. from reference allele
# -rmTriallelic 1e-5: remove high likelihood triallelic sites (p=value< .00001)
# -minInd 50: only keep sites with information (at least one read) from 50 individuals

#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J SAFAngsd
#SBATCH -o /home/ecalfee/hilo/slurm-log/SAFAngsd_%j_%A_%a.out
#SBATCH -t 36:00:00
#SBATCH --mem=60G
#SBATCH --array=1000,2000,3000,18-31,33-35,360-363,365-374

# pops >=1000 are the larger groupings of maize sympatric, mexicana sympatric and mexicana allopatric inds

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load angsd
module load angsd

# make directory to store output (if doesn't yet exist)
mkdir -p SFS/pass1/

echo "finding allele frequencies for pop$SLURM_ARRAY_TASK_ID"
# steps:
# (0) Start with filtered BAM files and reference genome
# (1) Use all sites to estimate site allele frequency
angsd -out SFS/pass1/pop$SLURM_ARRAY_TASK_ID \
-anc /group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa -fold 1 \
-bam pass1_bam_pops/pop$SLURM_ARRAY_TASK_ID.list \
-remove_bads 1 -minMapQ 30 \
-GL 1 -dosaf 1
echo "done with AFS"

# settings:
# -doSaf 1: Calculate the Site allele frequency likelihood based on individual genotype likelihoods assuming HWE
# -anc reference_genome with -fold 1 calculates a folded SFS (in lieue of having an ancestral state/genome)
# -remove_bads removes reads with flags like duplicates & -minMapQ 30: filter out sites with low mapping quality
# -bam list of bams to include
# -GL 1: use samtools genotype likelihood method 



#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J SAFByChr
#SBATCH -o /home/ecalfee/hilo/slurm-log/SAFByChr_%j_%A_%a.out
#SBATCH -t 64:00:00
#SBATCH --mem=25G
#SBATCH -n 3
#SBATCH --array=1-10

# POP variable needs to be set when calling the script
# e.g. sbatch --export=POP=18 calcSAFAngsd_byChr.sh
# pops >=1000 are the larger groupings of maize sympatric, mexicana sympatric and mexicana allopatric inds

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# rather than load module, use local module of angsd & print version
# angsd --version

# make directory to store output (if doesn't yet exist)
mkdir -p SFS/pass1/pop$POP

echo "finding allele frequencies for pop$POP, chr$SLURM_ARRAY_TASK_ID"
# steps:
# (0) Start with filtered BAM files and reference genome
# (1) Use all sites to estimate site allele frequency
angsd -out SFS/pass1/pop$POP/chr$SLURM_ARRAY_TASK_ID \
-anc /group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa -fold 1 \
-r $SLURM_ARRAY_TASK_ID \
-bam pass1_bam_pops/pop$POP.list \
-remove_bads 1 -minMapQ 30 \
-GL 1 -dosaf 1 \
-P 3
echo "done with SAF"

# settings:
# -r sets a region of the chromosome to do the SAF calculation -- here split by chromosomes
# -doSaf 1: Calculate the Site allele frequency likelihood based on individual genotype likelihoods assuming HWE
# -anc reference_genome with -fold 1 calculates a folded SFS (in lieue of having an ancestral state/genome)
# -remove_bads removes reads with flags like duplicates & -minMapQ 30: filter out sites with low mapping quality
# -bam list of bams to include
# -GL 1: use samtools genotype likelihood method
# -P sets number of angsd threads to use (but multithreading doesn't apply to IO read/write tasks)

#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J calcGL
#SBATCH -o /home/ecalfee/hilo/slurm-log/calcGL_%j_%A_%a.out
#SBATCH -t 5:00:00
#SBATCH --mem=15G
#SBATCH -n 4

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load angsd
module load angsd

# make directory to store output (if doesn't yet exist)
mkdir -p geno_lik/pass1/pruned_chunks/

# set variables
# pad the task id with leading zeros
printf -v TASK_ID "%03g" $SLURM_ARRAY_TASK_ID
# necessary becuase slurm array task ID doesn't hold leading zeros
POS_FILE=var_sites/pass1/pruned/positions_chunk$TASK_ID

# apply filtering with SAMtools & PICARD
echo "finding genotype likelihood using ANGSD on BAMS for hilo pruned SNPs chunk $SLURM_ARRAY_TASK_ID"
# steps:
# (0) Start with filtered BAM files, reference genome and list of positions in each chunk
# (1) For each chunk, index positions with ANGSD
angsd sites index $POS_FILE
# (2) make a list of chromosomes with positions in that chunk
cut -f1 $POS_FILE | uniq > $POS_FILE.chr
# (3) calculate genotype likelihoods using samtools algorithm and quality filters

angsd -out geno_lik/pass1/pruned_chunks/chunk_$TASK_ID \
-doMajorMinor 3 \
-sites $POS_FILE -rf $POS_FILE.chr \
-GL 1 -doGlf 2 \
-minMapQ 30 -minQ 20 \
-bam pass1_bam.all.list \
-remove_bads 1 \
-P 4

# settings:
# -sites and -rf specify which positions and chromosomes to calculate GL's for
# -remove_bads removes reads with flags like duplicates
# -doMajorMinor 3 takes major and minor allele from -sites file (for consistency)
# -bam list of bams to include (all newly sequenced allopatric mex. and sympatric mexicana & maize pops)
# -GL 1: use samtools genotype likelihood method
# -doGlf 2: output beagle likelihood file
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# -P 4 means use 4 threads/nodes for each angsd task

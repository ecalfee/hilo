#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J countReads
#SBATCH -o /home/ecalfee/hilo/slurm-log/countReadsFiltered_%A_%a.out
#SBATCH -t 24:00:00
#SBATCH --mem=100M
#SBATCH --array=1-40,44-79,81-200

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load samtools for quality read filtering
module load samtools

# make directory to store output (if doesn't yet exist)
mkdir -p countReads

# count raw reads after initial alignment
echo "flagstat original aln.bam file:" > countReads/hilo_$SLURM_ARRAY_TASK_ID.counts
samtools flagstat /group/jrigrp6/DanAlignments/HILO$SLURM_ARRAY_TASK_ID/aln.bam >> countReads/hilo_$SLURM_ARRAY_TASK_ID.counts

echo 'pre-filtering read count complete for hilo'$SLURM_ARRAY_TASK_ID

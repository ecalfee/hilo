#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J countReads
#SBATCH -o /home/ecalfee/hilo/slurm-log/countReadsFiltered_%j_%A_%a.out
#SBATCH -t 12:00:00
#SBATCH --mem=100M
#SBATCH --array=1-40,44-79,81-200

# to run (e.g. hilo1-hilo40): sbatch --array=1-40 countReadsFiltered.sh 

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

# append count only for reads passing mapping quality filtering
echo "nreads pass mapQ >= 30 filter: " $(samtools view -q 30 /group/jrigrp6/DanAlignments/HILO$SLURM_ARRAY_TASK_ID/aln.bam | wc -l) >> countReads/hilo_$SLURM_ARRAY_TASK_ID.counts

# append count for after de-duplication
echo "nreads after mapQ >= 30 and de-duplication filters: " $(samtools view -q 30 filtered_bam/hilo_$SLURM_ARRAY_TASK_ID.sort.dedup.bam | wc -l) >> countReads/hilo_$SLURM_ARRAY_TASK_ID.counts

echo 'pre- and post- filtering read count complete for hilo'$SLURM_ARRAY_TASK_ID

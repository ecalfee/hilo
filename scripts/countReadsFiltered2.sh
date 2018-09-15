#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J count2
#SBATCH -o /home/ecalfee/hilo/slurm-log/count2_%A_%a.out
#SBATCH -t 30:00:00
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
mkdir -p countReads2

# count raw reads after initial alignment
echo "total reads in original aln.bam file:" > countReads2/hilo_$SLURM_ARRAY_TASK_ID.raw.counts
samtools view -c /group/jrigrp6/DanAlignments/HILO$SLURM_ARRAY_TASK_ID/aln.bam >> countReads2/hilo_$SLURM_ARRAY_TASK_ID.raw.counts

# count mapped reads only
echo "mapped reads in original aln.bam file:" > countReads2/hilo_$SLURM_ARRAY_TASK_ID.mapped.counts
samtools view -c -F 4 /group/jrigrp6/DanAlignments/HILO$SLURM_ARRAY_TASK_ID/aln.bam >> countReads2/hilo_$SLURM_ARRAY_TASK_ID.mapped.counts

# count reads passing mapping quality filtering
echo "nreads pass mapQ >= 30 filter: " > countReads2/hilo_$SLURM_ARRAY_TASK_ID.mapQ30.counts
samtools view -c -q 30 /group/jrigrp6/DanAlignments/HILO$SLURM_ARRAY_TASK_ID/aln.bam >> countReads2/hilo_$SLURM_ARRAY_TASK_ID.mapQ30.counts

# count reads after de-duplication & map quality filtering
echo "nreads after mapQ >= 30 and de-duplication filters: " > countReads2/hilo_$SLURM_ARRAY_TASK_ID.mapQ30.dedup.counts
samtools view -c -q 30 filtered_bam/hilo_$SLURM_ARRAY_TASK_ID.sort.dedup.bam >> countReads2/hilo_$SLURM_ARRAY_TASK_ID.mapQ30.dedup.counts

echo 'read count complete for hilo'$SLURM_ARRAY_TASK_ID

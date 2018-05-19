#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo
#SBATCH -J hiloFiltBam
#SBATCH -o /home/ecalfee/hilo/slurm-log/out-%j_%a.txt
#SBATCH -e /home/ecalfee/hilo/slurm-log/error-%j_%a.txt
#SBATCH -t 0:05:00
#SBATCH --array=1-40

# load samtools for quality read filtering
module load samtools
# load picardtools for removing PCR duplicate reads
module load picardtools
# make directory to store output (if doesn't yet exist)
mkdir -p data/filtered_bam
# apply filtering with SAMtools
# I need to sort, index, remove duplicates, check that I am using samtools 1.6 if I want ot removedup that way
# or use picard tools
# and take out low mapping quality (< 30) and low adjusted base quality (<20) bases
samtools view -b /group/jrigrp6/DanAlignments/HILO$SLURM_ARRAY_TASK_ID/aln.bam
> data/filtered_bam/hilo_$SLURM_ARRAY_TASK_ID.bam
echo "filtering BAM for hilo $SLURM_ARRAY_TASK_ID"

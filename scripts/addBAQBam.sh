#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/filtered_bam
#SBATCH -J addBAQBam
#SBATCH -o /home/ecalfee/hilo/slurm-log/addBAQBam_%j_%A_%a.out
#SBATCH -t 24:00:00

# this script caps base quality scores with samtools adjusted scores (BAQ)
# which accounts for higher uncertainty in alignment around indels.
# the purpose is to use BAQ outside of samtools, e.g. when calling SNPs with ANGSD
# or counting reads that cover a specific base

# to run (e.g. hilo1-hilo40): sbatch addBAQBam.sh --array=1-40

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load samtools for quality read filtering
module load samtools

# calculate BAQ with SAMtools
echo "calculating BAQ for hilo $SLURM_ARRAY_TASK_ID"
# steps:
# (0) Start with sorted and deduplicated bam file: hilo_X.sort.dedup.bam

# (1) SAMTOOLS calculate BAQ > name.sort.dedup.baq.bam
samtools calmd -bAr hilo_$SLURM_ARRAY_TASK_ID.sort.dedup.bam \
-ref /group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa > hilo_$SLURM_ARRAY_TASK_ID.sort.dedup.baq.bam 

echo "done calculating BAQ"
# (2) SAMTOOLS index
samtools index hilo_$SLURM_ARRAY_TASK_ID.sort.dedup.baq.bam
echo "done indexing"

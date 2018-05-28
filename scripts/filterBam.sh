#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J filterBam
#SBATCH -o /home/ecalfee/hilo/slurm-log/filterBam_%j_%A_%a.out
#SBATCH -t 8:00:00
#SBATCH -x bigmem1

# jobs were failing due to lack of scratch memory so I exclude nodes with <30G
# available and I additionally exlude node bigmemm2 due to permissions issues

# to run (e.g. hilo1-hilo40): sbatch filterBam.sh --array=1-40

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# variable to keep track of any failed parts of the script
fails=false

# load samtools for quality read filtering
module load samtools
# load java to run picardtools
module load java
# load picardtools for removing PCR duplicate reads
module load picardtools # saves path to loaded versin in $PICARD variable
# make directory to store output (if doesn't yet exist)
mkdir -p filtered_bam || fails=true
# make a local ‘scratch’ directory for temporary files (@ end check that it’s empty)
mkdir -p /scratch/ecalfee/hilo_$SLURM_ARRAY_TASK_ID || fails=true   # temporary sort files will be written to local node

# apply filtering with SAMtools & PICARD
echo "filtering BAM for hilo $SLURM_ARRAY_TASK_ID"
# steps:
# (0) Start with raw bam file DIRECTORY/HILOX/aln.bam

# (1) SAMTOOLS sort reads by coordinate position > name.sort.bam
samtools sort -m 6G -T /scratch/ecalfee/hilo_$SLURM_ARRAY_TASK_ID \
-o /scratch/ecalfee/hilo_$SLURM_ARRAY_TASK_ID.sort.bam \
/group/jrigrp6/DanAlignments/HILO$SLURM_ARRAY_TASK_ID/aln.bam || fails=true

ls /scratch/ecalfee/
# (2) Picard MarkDuplicate marks and removes PCR duplicates (pipes directly to next step)
# note that –Xmx8G means anything over 8G memory will be written to the temporary directory TMP_DIR

# (3) SAMTOOLS removes low mapping quality reads (<30) > name.sort.dedup.bam
java -Xmx6g -jar $PICARD/picard.jar MarkDuplicates \
INPUT=/scratch/ecalfee/hilo_$SLURM_ARRAY_TASK_ID.sort.bam OUTPUT=/dev/stdout QUIET=true \
REMOVE_DUPLICATES=true TMP_DIR=/scratch/ecalfee/hilo_$SLURM_ARRAY_TASK_ID \
METRICS_FILE=filtered_bam/hilo_$SLURM_ARRAY_TASK_ID.metrics.txt | samtools view -b -q 30 \
-o filtered_bam/hilo_$SLURM_ARRAY_TASK_ID.sort.dedup.bam - || fails=true

# (4) remove intermediate file and temporary scratch directory
rm /scratch/ecalfee/hilo_$SLURM_ARRAY_TASK_ID.sort.bam
rm -r /scratch/ecalfee/hilo_$SLURM_ARRAY_TASK_ID # remove temporary scratch folder

# print confirmation that all parts ran without errors
if ! $fails; then
    echo 'NoFail - this script ran without errors on hilo $SLURM_ARRAY_TASK_ID'
fi

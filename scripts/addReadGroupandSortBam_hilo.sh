#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J RG+Sort
#SBATCH -o /home/ecalfee/hilo/slurm-log/addRGandSortBam_%A_%a.out
#SBATCH -t 2-00:00:00
#SBATCH --mem=16G
#SBATCH --array=1-200
#SBATCH --export=SEQ_LIBRARY=March2018


# this script pre-processes bams so that they are compatible with GATK:
# as the 1st step of bam file filtering, it adds read groups and sorts.
# The 2nd step is removing duplicates and calculating BAQ: dedupAndAddBAQ_hilo.sh

# array is the individual hilo id #
i=$SLURM_ARRAY_TASK_ID
DIR_IN="hilo_bam_mapped2v4_allContigs"
DIR_TMP="/scratch/ecalfee/sort_hilo_$i"
DIR_OUT="${DIR_IN}/results"

mkdir -p "${DIR_OUT}"
mkdir -p "${DIR_TMP}"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load samtools for quality read filtering
module load samtools
# load java to run picardtools
module load java
# load picardtools for removing PCR duplicate reads
module load picardtools # saves path to loaded versin in $PICARD variable

echo "adding Read Groups to BAM and sorting for hilo "$i
# (1) Add read group to BAM header using Picard
# (2) Sort reads by coordinate with samtools
java -jar $PICARD/picard.jar AddOrReplaceReadGroups \
INPUT="${DIR_IN}/hilo_$i.bam" \
OUTPUT=dev/stdout \
QUIET=true \
RGPL=illumina \
RGPU=HILO$i \
RGSM=HILO$i \
RGLB=$SEQ_LIBRARY | \
samtools sort -m 6G -T "${DIR_TMP}" \
-o "${DIR_OUT}/hilo_$i.sort.bam"

echo "all done!"

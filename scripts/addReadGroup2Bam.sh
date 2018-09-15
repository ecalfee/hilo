#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J addRG2Bam
#SBATCH -o /home/ecalfee/hilo/slurm-log/addRG2Bam_%A_%a.out
#SBATCH -t 2:00:00
#SBATCH -x bigmem1,bigmem7
#SBATCH --mem=16G
#SBATCH --array=1-200
#SBATCH --export=SEQ_LIBRARY=March2018

# array is the individual hilo id #
i=$SLURM_ARRAY_TASK_ID

# this script pre-processes bams so that they are compatible with GATK ASEReadCounter
# it adds read groups using picard and then re-indexes the new bams with samtools

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

echo "adding Read Groups to BAM for hilo "$i
# (1) Add read group to BAM header using Picard
java -jar $PICARD/picard.jar AddOrReplaceReadGroups I=filtered_bam/hilo_$i.sort.dedup.baq.bam \
O=filtered_bam/hilo_$i.sort.dedup.baq.rg.bam \
RGPL=illumina \
RGPU=HILO$i \
RGSM=HILO$i \
RGLB=$SEQ_LIBRARY


echo "indexing with SAMTOOLS for hilo "$i
# (2) index with SAMTOOLS
samtools index filtered_bam/hilo_$i.sort.dedup.baq.rg.bam

echo "all done!"

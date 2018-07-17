#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J reorderBam
#SBATCH -o /home/ecalfee/hilo/slurm-log/reorderBam_%j_%A_%a.out
#SBATCH -t 12:00:00
#SBATCH --mem=12G

# jobs will fail due to lack of scratch memory so I can exclude nodes with low available memory as needed, e.g. -x bigmem1

# to run (e.g. 4 ind's per lowland maize pop): sbatch reorderBam.sh --array=1-4,11-14,21-24,31-34

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# pad the task id with leading zeros
printf -v TASK_ID "%02g" $SLURM_ARRAY_TASK_ID
# necessary becuase slurm array task ID doesn't hold leading zeros

# load samtools for quality read filtering
module load samtools
# load java to run picardtools
module load java
# load picardtools for removing PCR duplicate reads
module load picardtools # saves path to loaded versin in $PICARD variable

# make a local ‘scratch’ directory for temporary files (@ end check that it’s empty)
mkdir -p /scratch/ecalfee/maizeLow_$TASK_ID   # temporary sort files will be written to local node

# reorder chromosomes to match APG4 reference
java -Xmx11g -jar $PICARD/picard.jar ReorderSam \
          INPUT=/alloMaize4pop_symlink_bam/maizeLow_$TASK_ID.bam \
          OUTPUT=/alloMaize4pop_symlink_bam/maizeLow_$TASK_ID.reordered.bam \
          TMP_DIR=/scratch/ecalfee/maizeLow_$TASK_ID \
          REFERENCE=/group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa \
          CREATE_INDEX=true
# (2) Picard ReorderSam re-orders a bam file to the chromosome ordering of a new (but same alignment) reference genome
# note that –Xmx8G means anything over 8G memory will be written to the temporary directory TMP_DIR
# CREATE_INDEX tells picard to make a new index file bam.bai for newly sorted bam


# (4) remove intermediate file and temporary scratch directory
rm -r /scratch/ecalfee/maizeLow_$TASK_ID # remove temporary scratch folder

# print confirmation that all parts ran without errors
echo "done re-ordering maizeLow "$TASK_ID" to APGv4"


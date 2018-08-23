#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J addBAQBamAllo
#SBATCH -o /home/ecalfee/hilo/slurm-log/addBAQBamAlloMaize_%j_%A_%a.out
#SBATCH -t 12:00:00
#SBATCH --mem=16G
#SBATCH --array=0-14

# this script caps base quality scores with samtools adjusted scores (BAQ)
# which accounts for higher uncertainty in alignment around indels.
# the purpose is to use BAQ outside of samtools, e.g. when calling SNPs with ANGSD
# or counting reads that cover a specific base

# it runs on the allopatric maize individuals used as a reference panel in the hilo study
# to run: sbatch addBAQBamAlloMaize.sh

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# variable to keep track of any failed parts of the script
LIST_OF_INDS=($(awk '{print $1}' landraces_fromLi/alloMaizeInclude.list)) # make array of individuals
IND=${LIST_OF_INDS[$SLURM_ARRAY_TASK_ID]} # get individual i based on SLURM_ARRAY_TASK_ID
INPUT_DIR="landraces_fromLi/remapped/" #redefine your input and output directory
OUTPUT_DIR="landraces_fromLi/remapped/"

# make temporary and output directories
mkdir -p $OUTPUT_DIR

# load samtools for quality read filtering
module load samtools

# calculate BAQ with SAMtools
echo "calculating BAQ for allopatric maize "$IND
# steps:
# (0) Start with sorted and deduplicated bam file: hilo_X.sort.dedup.bam

# (1) SAMTOOLS calculate BAQ > name.sort.dedup.baq.bam
samtools calmd -bArE --reference refMaize/AGPv4.fa \
$INPUT_DIR/$IND.sort.dedup.bam >> $OUTPUT_DIR/$IND.sort.dedup.bam
# -b bam output, # -Ar compute BAQ and cap base quality with BAQ 
# -E use extended BAQ calculation; this was a change to improve sensitivity and will be the new samtools default

echo "done calculating BAQ; now indexing allopatric maize "$IND
# (2) SAMTOOLS index
samtools index $OUTPUT_DIR/$IND.sort.dedup.bam
echo "done indexing allopatric maize "$IND

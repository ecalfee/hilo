#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J filterBamAllo
#SBATCH -o /home/ecalfee/hilo/slurm-log/filterBamAllo_%j_%A_%a.out
#SBATCH -t 24:00:00
#SBATCH --mem=16G
#SBATCH --array=0-14

# to run: sbatch filterBamAlloMaize.sh

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# variable to keep track of any failed parts of the script
LIST_OF_INDS=($(awk '{print $1}' landraces_fromLi/alloMaizeInclude.list)) # make array of individuals
IND=${LIST_OF_INDS[$SLURM_ARRAY_TASK_ID]} # get individual i based on SLURM_ARRAY_TASK_ID
TMPDIR="landraces_fromLi/tmp"$SLURM_ARRAY_TASK_ID #redefine your temporary directory here, you could make a new one in your work directory
INPUT_DIR="landraces_fromLi/remapped/" #redefine your input and output directory
OUTPUT_DIR="landraces_fromLi/remapped/"

# make temporary and output directories
mkdir -p $TMPDIR
mkdir -p $OUTPUT_DIR

# load samtools for quality read filtering
module load samtools
# load java to run picardtools
module load java
# load picardtools for removing PCR duplicate reads
module load picardtools # saves path to loaded versin in $PICARD variable

# apply filtering with SAMtools & PICARD
echo "done filtering BAM for allopatric maize "$IND
# steps:
# (0) Start with raw bam file

echo "sorting with SAMTOOLS "$IND
# (1) SAMTOOLS sort reads by coordinate position > name.sort.bam
samtools sort -m 14G -T $TMPDIR \
-o $OUTPUT_DIR/$IND.sort.bam \
$INPUT_DIR/$IND.bam

ls $TIMPDIR
ls $INPUT_DIR
# (2) Picard MarkDuplicate marks and removes PCR duplicates (pipes directly to next step)
# note that –Xmx8G means anything over 8G memory will be written to the temporary directory TMP_DIR
echo "removing duplicates with Picard and removing low mapping quality reads (mapQ < 30) with samtools "$IND
# (3) SAMTOOLS removes low mapping quality reads (<30) > name.sort.dedup.bam
java -Xmx6g -jar $PICARD/picard.jar MarkDuplicates \
INPUT=$OUTPUT_DIR/$IND.sort.bam OUTPUT=/dev/stdout QUIET=true \
REMOVE_DUPLICATES=true TMP_DIR=$TMPDIR \
METRICS_FILE=$OUTPUT_DIR/$IND.metrics.txt | samtools view -b -q 30 \
-o $OUTPUT_DIR/$IND.sort.dedup.bam -

# (4) remove intermediate file and temporary scratch directory
rm $OUTPUT_DIR/$IND.sort.bam
rm -r $TMPDIR # remove temporary scratch folder

# print confirmation that all parts ran without errors
echo "done filtering BAM for allopatric maize "$IND

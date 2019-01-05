#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J dedupBAQ
#SBATCH -o /home/ecalfee/hilo/slurm-log/dedupAndAddBAQ_%A_%a.out
#SBATCH -t 3-00:00:00

# to run (e.g. hilo1-hilo40): sbatch --array=1-40 dedupAndAddBAQ_hilo.sh

# this script takes in a sorted bam file, from addReadGroupandSortBam_hilo.sh
# and outputs a new bam with duplicates removed, BAQ scores calculated,
# and reads with low mapping quality (<30) filtered out. Also makes an index for new bam.

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

# command line arguments:
# note: all paths relative to bees/filtered_bams/
i=$SLURM_ARRAY_TASK_ID
DIR_MAIN="hilo_bam_mapped2v4_allContigs"
DIR_TMP="/scratch/ecalfee/dedup_hilo_$i" # for memory overflow
DIR_OUT="${DIR_MAIN}/results" # results directory
DIR_METRICS="${DIR_MAIN}/metrics" # metrics directory
BAM_SORTED="${DIR_OUT}/hilo_${i}.sort.bam" # full path to starting bam
BAM_OUT="${DIR_OUT}/hilo_${i}.sort.dedup.baq.bam"
REF="refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"

# make results directory (if necessary)
mkdir -p "${DIR_OUT}"
mkdir -p "${DIR_TMP}"
mkdir -p "${DIR_METRICS}"

echo "working directory: ${PWD}" # print current working directory
echo "picard path: ${PICARD}"

echo "marking duplicates with PICARD and calculating BAQ with SAMTOOLS"
java -Xmx6g -jar ${PICARD}/picard.jar MarkDuplicates \
INPUT="${BAM_SORTED}" OUTPUT=/dev/stdout QUIET=true \
REMOVE_DUPLICATES=true \
TMP_DIR="${DIR_TMP}" \
METRICS_FILE="${DIR_METRICS}/hilo_${i}.metrics.txt" | \
samtools calmd -SArE --reference ${REF} - | \
samtools view -bS -q 30 - > ${BAM_OUT}

echo "all done removing duplicates and calculating BAQ, now indexing!"
sleep 5s # because index needs to have a later timestamp
samtools index ${BAM_OUT}
echo "done indexing"

# options:
# -Xmx6g will not spill to tmp until 6G of memory are used
# samtools calmd calculates adjusted base quality score (BAQ)
# -S specifies that input is sam, not binary bam
# - sets stdin as input (stdout as output is default)
# -Ar compute BAQ and cap base quality with BAQ
# -E use extended BAQ calculation; this was a change to improve sensitivity and will be the new samtools default

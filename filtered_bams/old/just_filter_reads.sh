#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/filtered_bams
#SBATCH -J justFilt
#SBATCH -o /home/ecalfee/hilo/slurm-log/justFilter_%A_%a.out
#SBATCH -t 2-00:00:00
#SBATCH --mem=32G
#SBATCH -n 4
#SBATCH --array=0
#SBATCH --export=DIR_IN=../data/HILO_raw_reads,PREFIX_LIST=Jan2019

# START WITH ZERO FOR ARRAY INDEXES (!)

# this script takes in bam file with reads to APGv4 reference with all chromosome and contigs
# and sorted. It marks duplicates, calculates BAQ scores, and indexes the resulting BAM.
# to run: sbatch --array=1 --export=DIR_IN=../data/HILO_raw_reads/TEST,PREFIX_LIST=Jan2019 just_filter_reads.sh

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
# load bio for samtools and bwa
module load bio
# load java to run picardtools
module load java
# load picardtools for removing PCR duplicate reads
module load picardtools # saves path to loaded versin in $PICARD variable

echo "working directory: ${PWD}" # print current working directory
echo "picard path: ${PICARD}"

REF="../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
i=$SLURM_ARRAY_TASK_ID
ids=($(cat ${DIR_IN}/${PREFIX_LIST}_IDs.list))
ID=${ids["$i"]}
libraries=($(cat ${DIR_IN}/${PREFIX_LIST}_libraries.list))
LIBRARY=${libraries["$i"]}
lanes=($(cat ${DIR_IN}/${PREFIX_LIST}_lanes.list))
LANE=${lanes["$i"]}
DIR_OUT=results/${LANE} # results directory for final and intermediate bam files
DIR_METRICS=metrics/${LANE} # metrics directory
DIR_TMP="/scratch/ecalfee/${ID}" # for memory overflow

# make output directory
mkdir -p ${DIR_OUT}
mkdir -p ${DIR_METRICS}
mkdir -p ${DIR_TMP}

# note: this script is for reads already aligned to the reference genome & sorted

echo "marking duplicates with PICARD and calculating BAQ with SAMTOOLS"
java -Xmx6g -jar ${PICARD}/picard.jar MarkDuplicates \
INPUT="${DIR_OUT}/${ID}.sort.bam" OUTPUT=/dev/stdout QUIET=true \
REMOVE_DUPLICATES=true \
TMP_DIR="${DIR_TMP}" \
METRICS_FILE="${DIR_METRICS}/${ID}.metrics.txt" | \
samtools calmd -SArE --reference ${REF} - | \
samtools view -bS -q 30 - > "${DIR_OUT}/${ID}.sort.dedup.baq.bam"

echo "all done removing duplicates and calculating BAQ, now indexing!"
sleep 5s # because index needs to have a later timestamp
samtools index "${DIR_OUT}/${ID}.sort.dedup.baq.bam"
echo "done indexing"

echo "all done!"

# options

# picard markduplicates:
# -Xmx6g will not spill to tmp until 6G of memory are used

# samtools calmd calculates adjusted base quality score (BAQ):
# -S specifies that input is sam, not binary bam
# - sets stdin as input (stdout as output is default)
# -Ar compute BAQ and cap base quality with BAQ
# -E use extended BAQ calculation; this was a change to improve sensitivity and will be the new samtools default

#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/analyses/filtered_bams
#SBATCH -J mapFilt
#SBATCH -o /home/ecalfee/hilo/slurm-log/mapAndFilter_%A_%a.out
#SBATCH -t 6-00:00:00
#SBATCH --mem=32G
#SBATCH -n 16
#SBATCH --array=0-6
#SBATCH --export=DIR_IN=../../data/HILO_raw_reads,PREFIX_LIST=Jan2019

# START WITH ZERO FOR ARRAY INDEXES (!)

# this script takes in fastq files with raw paired reads
# and aligns reads to APGv4 reference with all chromosome and contigs
# sample ID, fastq file 1, fastq file 2, output bam directory
# to run: sbatch --array=1-$(wc -l prefix_ids_and_lanes_IDs.list) map_and_filter_reads.sh directory/with/fastq_files prefix_ids_and_lanes

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

REF="refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
i=$SLURM_ARRAY_TASK_ID
ids=($(cat ${PREFIX_LIST}_IDs.list))
ID=${ids["$i"]}
libraries=($(cat ${PREFIX_LIST}_libraries.list))
LIBRARY=${libraries["$i"]}
lanes=($(cat ${PREFIX_LIST}_lanes.list))
LANE=${lanes["$i"]}
FASTQ1=${DIR_IN}/${LANE}/${ID}_1.fq.gz
FASTQ2=${DIR_IN}/${LANE}/${ID}_2.fq.gz
DIR_OUT=results/${LANE} # results directory for final and intermediate bam files
DIR_METRICS=metrics/${LANE} # metrics directory
DIR_TMP="/scratch/ecalfee/${ID}" # for memory overflow

# make output directory
mkdir -p ${DIR_OUT}
mkdir -p ${DIR_METRICS}
mkdir -p ${DIR_TMP}

# align fastq to maize reference AGPv4 official release using bwa mem
# note: reference genome needs to already be indexed, e.g. bwa index ref.fa
echo "running bwa mem for sample ${ID} lane ${LANE} library ${LIBRARY}"
bwa mem -t 16 -v 3 \
-R "@RG\tID:${LANE}\tSM:${ID}\tPL:ILLUMINA\tLB:${LIBRARY}\tPU:${LANE}.${ID}"
"${REF}" "${FASTQ1}" "${FASTQ2}"  | \
samtools view -bS -o "${DIR_OUT}/${ID}.bam" -

echo "sorting BAM file"
samtools sort -m 6G -@ 4 -T "${DIR_TMP}" \
-o "${DIR_OUT}/${ID}.sort.bam" "${DIR_OUT}/${ID}.bam"

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
# bwa mem:
# -t k specifies k threads
# -v 3 is for extra verbose (comments in addition to error messages)
# 2 input fastq files (for paired reads)
# samtools sort:
# - means read from standard in
# -bS means SAM in and BAM out
# -@ specifies number of threads and -m amount of memory per thread
# -T specifies temporary directory to write files if exceeding memory limit

# picard markduplicates:
# -Xmx6g will not spill to tmp until 6G of memory are used

# samtools calmd calculates adjusted base quality score (BAQ):
# -S specifies that input is sam, not binary bam
# - sets stdin as input (stdout as output is default)
# -Ar compute BAQ and cap base quality with BAQ
# -E use extended BAQ calculation; this was a change to improve sensitivity and will be the new samtools default

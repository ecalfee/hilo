#!/bin/bash
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/functional_analysis
#SBATCH -J 3fq2AB10
#SBATCH -o /home/ecalfee/hilo/slurm-log/3fastq2AB10_%A_%a.out
#SBATCH -t 6:00:00
#SBATCH --mem=32G
#SBATCH -n 4

# this script creates a bam mapped to ab10 and competitive regions starting with 3 fastq files --
# PREFIX_1.fq.gz PREFIX_2.fq.gz PREFIX_3.fg.gz
# where _1 and _2 are the paired reads and _3 is the set of singletons
# note: fastq files must be name (query) sorted, so read pairs can be identified
# after merging aligned bams, the script deduplicates, adds BAQ, and indexes final bam

# to run: sbatch -p bigmemm --array=0-14 --export=ID_DIR=../data/landraces_fromLi,DIR_IN=../data/landraces_fromLi/original/fastq,PREFIX_LIST=alloMAIZE threefastq2ab10.sh

# load bio for samtools and bwa
module load bio
# load java to run picardtools
module load java
# load picardtools for removing PCR duplicate reads
module load picardtools # saves path to loaded versin in $PICARD variable
echo "working directory: ${PWD}" # print current working directory
echo "picard path: ${PICARD}"

REF="../data/ab10/kinesins.fa"
i=$SLURM_ARRAY_TASK_ID
ids=($(cat ${ID_DIR}/${PREFIX_LIST}_IDs.list))
ID=${ids["$i"]}
libraries=($(cat ${ID_DIR}/${PREFIX_LIST}_libraries.list))
LIBRARY=${libraries["$i"]}
lanes=($(cat ${ID_DIR}/${PREFIX_LIST}_lanes.list))
LANE=${lanes["$i"]}
FASTQ1=${DIR_IN}/${ID}_1.fq.gz
FASTQ2=${DIR_IN}/${ID}_2.fq.gz
FASTQ3=${DIR_IN}/${ID}_3.fq.gz
DIR_OUT=results/bams_ab10/${LANE} # results directory for final and intermediate bam files
DIR_METRICS=metrics_ab10/${LANE} # metrics directory
DIR_TMP="/scratch/ecalfee/${ID}" # for memory overflow


# make output directory
mkdir -p ${DIR_OUT}
mkdir -p ${DIR_METRICS}
mkdir -p ${DIR_TMP}

# align fastq to ab10 and competitive sequences
# note: reference genome needs to already be indexed, e.g. bwa index ref.fa
echo "running bwa mem for sample ${ID} lane ${LANE} library ${LIBRARY} paired reads"
bwa mem -t 16 -v 3 \
-R "@RG\tID:${LANE}\tSM:${ID}\tPL:ILLUMINA\tLB:${LIBRARY}\tPU:${ID}.${LANE}" \
"${REF}" "${FASTQ1}" "${FASTQ2}"  | \
samtools view -Shu -q 1 - | \
samtools sort -m 6G -@ 4 -T "${DIR_TMP}" - > "${DIR_OUT}/${ID}_12.ab10.sort.bam"
# -q 1 filter only keeps reads that map at all

echo "running bwa mem for sample ${ID} lane ${LANE} library ${LIBRARY} UNpaired reads"
bwa mem -t 16 -v 3 \
-R "@RG\tID:${LANE}\tSM:${ID}\tPL:ILLUMINA\tLB:${LIBRARY}\tPU:${ID}.${LANE}" \
"${REF}" "${FASTQ3}"  | \
samtools view -Shu -q 1 - | \
samtools sort -m 6G -@ 4 -T "${DIR_TMP}" - > "${DIR_OUT}/${ID}_3.ab10.sort.bam"
# -q 1 filter only keeps reads that map at all

echo "merging paired and unpaired reads into final BAM"
samtools merge "${DIR_OUT}/${ID}.ab10.sort.bam" "${DIR_OUT}/${ID}_12.ab10.sort.bam" "${DIR_OUT}/${ID}_3.ab10.sort.bam"

echo "marking duplicates with PICARD and calculating BAQ with SAMTOOLS"
java -Xmx6g -jar ${PICARD}/picard.jar MarkDuplicates \
INPUT="${DIR_OUT}/${ID}.ab10.sort.bam" OUTPUT=/dev/stdout QUIET=true \
REMOVE_DUPLICATES=true \
TMP_DIR="${DIR_TMP}" \
METRICS_FILE="${DIR_METRICS}/${ID}.ab10.metrics.txt" | \
samtools calmd -SArE --reference ${REF} - | \
samtools view -bS -q 30 - > "${DIR_OUT}/${ID}.ab10.sort.dedup.baq.bam"

# note: I do not delete intermediate bam (should be quite small)

echo "all done removing duplicates and calculating BAQ, now indexing!"
sleep 5s # because index needs to have a later timestamp
samtools index "${DIR_OUT}/${ID}.ab10.sort.dedup.baq.bam"
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


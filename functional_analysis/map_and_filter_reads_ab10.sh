#!/bin/bash
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/functional_analysis
#SBATCH -J mapAB10
#SBATCH -o /home/ecalfee/hilo/slurm-log/mapAB10_%A_%a.out
#SBATCH -t 6:00:00
#SBATCH --mem=32G
#SBATCH -n 4
#SBATCH --array=0
#SBATCH --export=DIR_IN=../data/HILO_raw_reads,PREFIX_LIST=Jan2019

# START WITH ZERO FOR ARRAY INDEXES (!)

# this script takes in fastq files with raw paired reads
# and aligns reads to ab10 fasta (with competitive making to other similar regions)
# to run: sbatch --array=1 --export=DIR_IN=../data/HILO_raw_reads/TEST,PREFIX_LIST=Jan2019 map_and_filter_reads_ab10.sh

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

REF="../data/incompatibility_loci/kinesins.fa"
i=$SLURM_ARRAY_TASK_ID
ids=($(cat ${DIR_IN}/${PREFIX_LIST}_IDs.list))
ID=${ids["$i"]}
libraries=($(cat ${DIR_IN}/${PREFIX_LIST}_libraries.list))
LIBRARY=${libraries["$i"]}
lanes=($(cat ${DIR_IN}/${PREFIX_LIST}_lanes.list))
LANE=${lanes["$i"]}
FASTQ1=${DIR_IN}/${LANE}/${ID}_1.fq.gz
FASTQ2=${DIR_IN}/${LANE}/${ID}_2.fq.gz
DIR_OUT=results/bams_ab10/${LANE} # results directory for final and intermediate bam files
DIR_METRICS=metrics_ab10/${LANE} # metrics directory
DIR_TMP="/scratch/ecalfee/${ID}" # for memory overflow

# make output directory
mkdir -p ${DIR_OUT}
mkdir -p ${DIR_METRICS}
mkdir -p ${DIR_TMP}

# align fastq to ab10 and competitive sequences
# note: reference genome needs to already be indexed, e.g. bwa index ref.fa
echo "running bwa mem for sample ${ID} lane ${LANE} library ${LIBRARY} and sorting BAM file"
bwa mem -t 16 -v 3 \
-R "@RG\tID:${LANE}\tSM:${ID}\tPL:ILLUMINA\tLB:${LIBRARY}\tPU:${ID}.${LANE}" \
"${REF}" "${FASTQ1}" "${FASTQ2}"  | \
samtools view -Shu -q 0 - | \
samtools sort -m 6G -@ 4 -T "${DIR_TMP}" - > "${DIR_OUT}/${ID}.ab10.sort.bam"

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
# -Ar compute BAQ and cap base quality with BAQ
# -E use extended BAQ calculation; this was a change to improve sensitivity and will be the new samtools default

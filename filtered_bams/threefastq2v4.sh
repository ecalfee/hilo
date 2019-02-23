#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/filtered_bams
#SBATCH -J 3fq2v4
#SBATCH -o /home/ecalfee/hilo/slurm-log/threefastq2v4_%A_%a.out
#SBATCH -t 10-00:00:00
#SBATCH --mem=64G
#SBATCH -n 16

# this script creates a bam mapped to v4 ref genome starting with 3 fastq files --
# PREFIX_1.fq.gz PREFIX_2.fq.gz PREFIX_3.fg.gz
# where _1 and _2 are the paired reads and _3 is the set of singletons
# note: fastq files must be name (query) sorted, so read pairs can be identified
# after merging aligned bams, the script deduplicates, adds BAQ, and indexes final bam

# to run: sbatch --export=ID=RIMMA0625,DIR_IN=/group/jrigrp4/landraces_v4_fromLi2017/original/fastq,LANE=alloMAIZE threefastq2v4.sh

DIR_TMP="/scratch/ecalfee/${ID}" # for memory overflow
REF="../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
DIR_OUT=results/${LANE} # results directory for final and intermediate bam files
DIR_METRICS=metrics/${LANE} # metrics directory

# so these .fq files can then be remapped to the APGv4 reference
mkdir -p "$DIR_OUT"
mkdir -p "$DIR_TMP"
mkdir -p "$DIR_METRICS"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load bio # loads samtools and bedtools
module load picardtools # saves path to loaded versin in $PICARD variable

echo "mapping paired reads with bwa for sample $ID and sorting with samtools"
# align fastq to maize reference AGPv4 official release using bwa mem
# note: reference genome needs to already be indexed, e.g. bwa index ref.fa
# IF bams already exist, they will be used rather than re-mapping (!)
if [ ! -f "${DIR_OUT}/${ID}_12.sort.bam" ]; then
	bwa mem -t 16 -v 3 \
	-R "@RG\tID:paired\tSM:${ID}\tPL:ILLUMINA\tLB:paired\tPU:${ID}" \
	"${REF}" "$DIR_IN"/"$ID"_1.fq.gz "$DIR_IN"/"$ID"_2.fq.gz | \
	samtools view -Shu - | \
	samtools sort -m 6G -@ 4 -T "${DIR_TMP}" - > "${DIR_OUT}/${ID}_12.sort.bam"
else
	echo "file ${DIR_OUT}/${ID}_12.sort.bam already exists and will be used"
fi

echo "mapping unpaired reads with bwa for sample $ID and sorting with samtools"

if [ ! -f "${DIR_OUT}/${ID}_3.sort.bam" ]; then
	bwa mem -t 16 -v 3 \
	-R "@RG\tID:unpaired\tSM:${ID}\tPL:ILLUMINA\tLB:unpaired\tPU:${ID}" \
	"${REF}" "$DIR_IN"/"$ID"_3.fq.gz  | \
	samtools view -Shu - | \
	samtools sort -m 6G -@ 4 -T "${DIR_TMP}" - > "${DIR_OUT}/${ID}_3.sort.bam"
else
	echo "file ${DIR_OUT}/${ID}_3.sort.bam already exists and will be used"
fi

echo "merging paired and unpaired reads into final BAM"
samtools merge "${DIR_OUT}/${ID}.sort.bam" "${DIR_OUT}/${ID}_12.sort.bam" "${DIR_OUT}/${ID}_3.sort.bam"

echo "marking duplicates with PICARD and calculating BAQ with SAMTOOLS"
java -Xmx6g -jar ${PICARD}/picard.jar MarkDuplicates \
INPUT="${DIR_OUT}/${ID}.sort.bam" OUTPUT=/dev/stdout QUIET=true \
REMOVE_DUPLICATES=true \
TMP_DIR="${DIR_TMP}" \
METRICS_FILE="${DIR_METRICS}/${ID}.metrics.txt" | \
samtools calmd -SArE --reference ${REF} - | \
samtools view -bS -q 30 - > "${DIR_OUT}/${ID}.sort.dedup.baq.bam"

echo "removing intermediate bam"
rm "${DIR_OUT}/${ID}.sort.bam"

echo "all done removing duplicates and calculating BAQ, now indexing!"
sleep 5s # because index needs to have a later timestamp
samtools index "${DIR_OUT}/${ID}.sort.dedup.baq.bam"
echo "done indexing"

echo "all done!"

# options
# samtools view
# -S is for sam input, -u for output uncompressed bam, -h is for keeping the header
# samtools sort
# -m is memory per thread, -@ is number of threads, and -T for temporary directory

# picard markduplicates:
# -Xmx6g will not spill to tmp until 6G of memory are used

# samtools calmd calculates adjusted base quality score (BAQ):
# -S specifies that input is sam, not binary bam
# - sets stdin as input (stdout as output is default)
# -Ar compute BAQ and cap base quality with BAQ
# -E use extended BAQ calculation; this was a change to improve sensitivity and will be the new samtools default



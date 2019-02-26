#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/filtered_bams
#SBATCH -J bam2fq
#SBATCH -o /home/ecalfee/hilo/slurm-log/bam2fastqMAIZE4LOW_%A_%a.out
#SBATCH -t 10-00:00:00
#SBATCH --mem=48G
#SBATCH -n 16
#SBATCH --array=1-4,11-14,21-24,31-34

# this script takes in a bam file and creates 3 fastq files --
# PREFIX_1.fq.gz PREFIX_2.fq.gz PREFIX_3.fg.gz
# where _1 and _2 are the paired reads and _3 is the set of singletons
# before creating the fastq files, reads must be sorted by query name so pairs are in the same order
# then the script creates a bam mapped to v4 ref genome starting with 3 fastq files --
# PREFIX_1.fq.gz PREFIX_2.fq.gz PREFIX_3.fg.gz
# where _1 and _2 are the paired reads and _3 is the set of singletons
# note: fastq files must be name (query) sorted, so read pairs can be identified
# after merging aligned bams, the script deduplicates, adds BAQ, and indexes final bam


# to run: sbatch --array=1-4,11-14,21-24,31-34 --export=DIR_IN=../data/alloMaize4Low,LANE=MAIZE4LOW bam2fastq2v4_MAIZE4LOW.sh

ID="MAIZE4LOW$SLURM_ARRAY_TASK_ID"
DIR_TMP="/scratch/ecalfee/${ID}" # for memory overflow
REF="../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
DIR_OUT=results/${LANE} # results directory for final and intermediate bam files
DIR_METRICS=metrics/${LANE} # metrics directory

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load bio # loads samtools and bedtools
module load picardtools # saves path to loaded versin in $PICARD variable

# make directories
mkdir -p "$DIR_IN"/fastq
mkdir -p "$DIR_OUT"
mkdir -p "$DIR_TMP"
mkdir -p "$DIR_METRICS"


echo "sorting BAM by name then output fastq"
# use samtools to sort input bam file by name (so paired reads are together)
# then use samtools to output fastq files
samtools sort -n "$DIR_IN"/"$BAM" | \
samtools fastq -c 6 -1 "$DIR_IN"/fastq/"$ID"_1.fq.gz -2 "$DIR_IN"/fastq/"$ID"_2.fq.gz -s "$DIR_IN"/fastq/"$ID"_3.fq.gz -

# align fastq to maize reference AGPv4 official release using bwa mem
echo "mapping paired reads with bwa for sample $ID and sorting with samtools"

bwa mem -t 16 -v 3 \
-R "@RG\tID:paired\tSM:${ID}\tPL:ILLUMINA\tLB:paired\tPU:${ID}" \
"${REF}" "$DIR_IN"/"$ID"_1.fq.gz "$DIR_IN"/"$ID"_2.fq.gz | \
samtools view -Shu - | \
samtools sort -m 6G -@ 6 -T "${DIR_TMP}" - > "${DIR_OUT}/${ID}_12.sort.bam"

echo "mapping unpaired reads with bwa for sample $ID and sorting with samtools"

bwa mem -t 16 -v 3 \
-R "@RG\tID:unpaired\tSM:${ID}\tPL:ILLUMINA\tLB:unpaired\tPU:${ID}" \
"${REF}" "$DIR_IN"/"$ID"_3.fq.gz  | \
samtools view -Shu - | \
samtools sort -m 6G -@ 6 -T "${DIR_TMP}" - > "${DIR_OUT}/${ID}_3.sort.bam"

echo "merging paired and unpaired reads into final BAM"
samtools merge "${DIR_OUT}/${ID}.sort.bam" "${DIR_OUT}/${ID}_12.sort.bam" "${DIR_OUT}/${ID}_3.sort.bam"

echo "removing temporary bams _12 and _3"
rm "${DIR_OUT}/${ID}_12.sort.bam" "${DIR_OUT}/${ID}_3.sort.bam"

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

# samtools sort
# -n is for sorting by read name

# samtools fastq
# -s is for singletons, -1 for read 1 and -2 for read 2
# -c is for compression, level 6 is the default for gzip

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


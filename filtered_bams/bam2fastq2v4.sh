#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/filtered_bams
#SBATCH -J bam2fq
#SBATCH -o /home/ecalfee/hilo/slurm-log/bam2fastq_%A_%a.out
#SBATCH -t 8-00:00:00
#SBATCH --mem=32G
#SBATCH -n 16

# this script takes in a bam file and creates 3 fastq files --
# PREFIX_1.fq.gz PREFIX_2.fq.gz PREFIX_3.fg.gz
# where _1 and _2 are the paired reads and _3 is the set of singletons
# before creating the fastq files, reads must be sorted by query name so pairs are in the same order

# to run: sbatch --export=ID=RIMMA0625,DIR_IN=../data/landraces_fromLi/original,BAM=RIMMA0625_Andean.IndelRealigned.bam bam2fastq2v4.sh

DIR_TMP="/scratch/ecalfee/${ID}" # for memory overflow
REF="../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"


# so these .fq files can then be remapped to the APGv4 reference
mkdir -p "$DIR_IN"/fastq
mkdir -p "$DIR_IN"/remapped
mkdir -p "$DIR_TMP"


# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load bio # loads samtools and bedtools

echo "sorting BAM by name then output fastq"
# use samtools to sort input bam file by name (so paired reads are together)
# then use samtools to output fastq files
samtools sort -n "$BAM" | \
samtools fastq -c 6 -1 "$DIR_IN"/fastq/"$ID"_1.fq.gz -2 "$DIR_IN"/fastq/"$ID"_2.fq.gz -s "$DIR_IN"/fastq/"$ID"_3.fq.gz -


echo "mapping paired reads with bwa for sample $ID"
# align fastq to maize reference AGPv4 official release using bwa mem
# note: reference genome needs to already be indexed, e.g. bwa index ref.fa
bwa mem -t 16 -v 3 \
-R "@RG\tID:paired\tSM:${ID}\tPL:ILLUMINA\tLB:paired\tPU:${ID}" \
"${REF}" "$DIR_IN"/fastq/"$ID"_1.fq.gz "$DIR_IN"/fastq/"$ID"_2.fq.gz  | \
samtools view -bS -o "${DIR_IN}/remapped/${ID}_12.bam" -

echo "mapping unpaired reads with bwa for sample $ID"
bwa mem -t 16 -v 3 \
-R "@RG\tID:unpaired\tSM:${ID}\tPL:ILLUMINA\tLB:unpaired\tPU:${ID}" \
"${REF}" "$DIR_IN"/fastq/"$ID"_3.fq.gz  | \
samtools view -bS -o "${DIR_IN}/remapped/${ID}_3.bam" -

echo "sorting BAM - paired reads"
samtools sort -m 6G -@ 4 -T "${DIR_TMP}" \
-o "${DIR_IN}/remapped/${ID}_12.sort.bam" "${DIR_IN}/remapped/${ID}_12.bam"

echo "sorting BAM - unpaired reads"
samtools sort -m 6G -@ 4 -T "${DIR_TMP}" \
-o "${DIR_IN}/remapped/${ID}_3.sort.bam" "${DIR_IN}/remapped/${ID}_3.bam"

echo "merging paired and unpaired reads into final BAM"
samtools merge "${DIR_IN}/remapped/${ID}.sort.bam" "${DIR_IN}/remapped/${ID}_12.sort.bam" "${DIR_IN}/remapped/${ID}_3.sort.bam"

echo "all done!"



# options

# samtools sort
# -n is for sorting by read name

# samtools fastq
# -s is for singletons, -1 for read 1 and -2 for read 2
# -c is for compression, level 6 is the default for gzip

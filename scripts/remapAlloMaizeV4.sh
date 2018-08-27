#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J remapV4
#SBATCH -o /home/ecalfee/hilo/slurm-log/remapV4_%j_%A_%a.out
#SBATCH -t 192:00:00
#SBATCH --mem=24G
#SBATCH -n 16
#SBATCH --array=0-14

# this script takes bam files with reads aligned to APGv4 reference with
# many unassigned contigs -> convertes to fastq -> then realigns to the APGv4 reference
# with only autosomes 1-10 and mt/pt

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load samtools
module load bwa

LIST_OF_INDS=($(awk '{print $1}' landraces_fromLi/alloMaizeInclude.list)) # make array of individuals
IND=${LIST_OF_INDS[$SLURM_ARRAY_TASK_ID]} # get individual i based on SLURM_ARRAY_TASK_ID
TMPDIR="landraces_fromLi/tmp" #redefine your temporary directory here, you could make a new one in your work directory
INPUT_DIR="landraces_fromLi/original/" #redefine your input and output directory
OUTPUT_DIR="landraces_fromLi/remapped/"

# make temporary and output directories
mkdir -p $TMPDIR
mkdir -p $OUTPUT_DIR

echo "BAM -> fastq for landrace "$IND
# use samtools to convert bam to fastq file only if fastq does not already exist
if [ ! -e "${TMPDIR}/${IND}.fq" ]
then
    echo "fastq does not exist; making .fq file"
    samtools fastq -0 /dev/null -F 0x900 ${INPUT_DIR}/${IND}.bam > ${TMPDIR}/$IND.fq
else
    echo "fastq already exists: ${TMPDIR}/${IND}.fq"
fi


echo "running bwa mem for landrace "$IND
# realign fastq using bwa to maize reference AGPv4 with autosomes 1-10 and pt/mt, but no extra contigs
bwa mem -t 16 -p refMaize/AGPv4.fa ${TMPDIR}/${IND}.fq  > ${TMPDIR}/${IND}.sam
# -t k specifies k threads
# -p specifies that the paired end reads in .fq are interwoven:
# read1_pair1, read1_pair2, read2_pair1, read2_pair2 etc.

echo "SAM -> BAM for landrace "$IND
samtools view -bS -o ${OUTPUT_DIR}/${IND}.bam ${TMPDIR}/${IND}.sam

#echo "deleting intermediate SAM & fastq files for landrace "$IND
#you might want to clear the sam files and fastq files, those take a lot of space.
#rm ${TMPDIR}/${IND}.sam
#rm ${TMPDIR}/${IND}.fq # will delete after I see the script has run properly

echo "all done for landrace "$IND

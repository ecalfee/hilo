#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J geneOver10kb
#SBATCH -o /home/ecalfee/hilo/slurm-log/geneOverlapWindows10KB_%A_%a.out
#SBATCH -t 2:00:00
#SBATCH --mem=8G

# script to divide genome into fixed non-overlapping 10kb windows 
# and calculate gene density within each window and recombination rates by window

CDS_FILE="refMaize/geneAnnotations/CDS_merged.bed"
GENOME_FILE="refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"
WINDOW_SIZE_BP=10000
DIR_OUT="refMaize/windows_10kb/"
WINDOWS_FILE=${DIR_OUT}"/whole_genome.bed"
GENE_OVERLAP_FILE=${DIR_OUT}"/gene_overlap.bed"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load software needed
module load bio
module load R

# make output directory
mkdir -p ${DIR_OUT}

# first divide genome into 10kb windows: whole_genome_10kb_windows.bed
echo "dividing genome into 10kb windows"
bedtools makewindows -g ${GENOME_FILE} -w ${WINDOW_SIZE_BP} -i winnum > ${WINDOWS_FILE}

# then calculate overlap:
echo "calculating gene overlap"
bedtools coverage -sorted -a ${WINDOWS_FILE} -b ${CDS_FILE} > ${GENE_OVERLAP_FILE}

# then calculate recombination rates for windows:
echo "calculate recombination rate"
Rscript ../scripts/bed_wind_2_recomb_rate.R ${WINDOWS_FILE}

echo 'all done!'

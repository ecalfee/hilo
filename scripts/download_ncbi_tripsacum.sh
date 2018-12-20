#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/
#SBATCH -J downTRIP
#SBATCH -o /home/ecalfee/hilo/slurm-log/downTRIP_%A_%a.out
#SBATCH -t 240:00:00
#SBATCH --mem=8G
#SBATCH -n 1
#SBATCH --array=0-1

# this script downloads raw reads from NCBI
# for high-coverage Tripsacum (outgroup) genome
# using 2 libraries: one mate-pair and one paired-end reads
# NCBI information: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN09910625

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load bio3

i=$SLURM_ARRAY_TASK_ID
LIST_OF_SRAs=($(cat tripsacum/SRR_Acc_List.txt))
SRA_ID=${LIST_OF_SRAs[$i]}
OUTPUT_DIR="tripsacum" # directory to store .fq files downloaded

# output directories
mkdir -p ${OUTPUT_DIR}

echo "Downloading data for SRA accession "${SRA_ID}
prefetch -v -X 90G -O ${OUTPUT_DIR} ${SRA_ID}

echo "all done for SRA accession "${SRA_ID}
# options:
# -v for verbose output
# -X is maximum file size to downloads in Gb -- these are big files!
# -O specifies output directory to save downloaded file

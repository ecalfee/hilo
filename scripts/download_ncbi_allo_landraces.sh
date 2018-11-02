#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/
#SBATCH -J downLANDRACE
#SBATCH -o /home/ecalfee/hilo/slurm-log/downLANDRACE_%A_%a.out
#SBATCH -t 12:00:00
#SBATCH --mem=16G
#SBATCH -n 4
#SBATCH --array=1-42

# array starts at 1 because skips first line (header) in SraRunTable.txt (see below)

# this script downloads raw reads from NCBI
# for allopatric maize landraces from
# Wang et al. Genome Biology (2017) 18:215 DOI 10.1186/s13059-017-1346-4

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load bio3

LIST_OF_SAMPLEs=($(awk '{print $12}' landraces_fromLi/SraRunTable.txt))
LIST_OF_SRAs=($(awk '{print $10}' landraces_fromLi/SraRunTable.txt))
SAMPLE_ID=${LIST_OF_SAMPLEs[$SLURM_ARRAY_TASK_ID]} # get individual i based on SLURM_ARRAY_TASK_ID
SRA_ID=${LIST_OF_SRAs[$SLURM_ARRAY_TASK_ID]}
OUTPUT_DIR="landraces_fromLi/ncbi/"${SAMPLE_ID} # directory to store .fq files downloaded

# output directories
mkdir -p ${OUTPUT_DIR}/tmp

echo "Downloading data for "${SAMPLE_ID}
echo "SRA accession "${SRA_ID}
parallel-fastq-dump --sra-id ${SRA_ID} \
--threads 4 --outdir ${OUTPUT_DIR} \
--gzip --clip --split-files --skip-technical \
--tmpdir ${OUTPUT_DIR}/tmp

echo "all done for landrace "${SAMPLE_ID}"; SRA accession "${SRA_ID}

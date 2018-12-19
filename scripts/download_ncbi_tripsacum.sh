#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/
#SBATCH -J downTRIP
#SBATCH -o /home/ecalfee/hilo/slurm-log/downTRIP_%A_%a.out
#SBATCH -t 120:00:00
#SBATCH --mem=4G
#SBATCH -n 4
#SBATCH --array=0-1

# this script downloads raw reads from NCBI
# for high-coverage Tripsacum (outgroup) genome
# using 2 libraries: one mate-pair and one paired-end reads

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
SCRATCH_DIR="/scratch/ecalfee/tripsacum_"$i
TMP_DIR=${SCRATCH_DIR}"/tmp"

# output directories
mkdir -p ${OUTPUT_DIR}
# scratch directories
mkdir -p ${SCRATCH_DIR}
mkdir -p ${TMP_DIR}

echo "Downloading data for SRA accession "${SRA_ID}
parallel-fastq-dump --sra-id ${SRA_ID} \
--threads 4 --outdir ${SCRATCH_DIR} \
--gzip --clip --split-files --skip-technical \
--tmpdir ${TMP_DIR}

# copy results back over
echo "all done downloading from NCBI! Now transfering files from local to home directory"
rsync -avh --remove-source-files ${DIR_SCRATCH}/ ${DIR_OUT}/ # copies all contents of output directory over to the appropriate home directory & cleans up scratch dir.
echo "results copied to home output directory: "${DIR_OUT}

echo "all done for SRA accession "${SRA_ID}

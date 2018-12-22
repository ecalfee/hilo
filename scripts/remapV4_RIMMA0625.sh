#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/
#SBATCH -J remapV4_RIMMA0625
#SBATCH -o /home/ecalfee/hilo/slurm-log/remapV4_RIMMA0625_%A_%a.out
#SBATCH -t 192:00:00
#SBATCH --mem=32G
#SBATCH -n 16

BAM_IN="landraces_fromLi/original/RIMMA0625_Andean.IndelRealigned"
ID="RIMMA0625"
DIR_OUT="landraces_fromLi/original"

# this script takes a bam with reads aligned to v3 -> fastq -> maps with bwa to v4 genome

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load samtools
module load bwa

# make output directory
mkdir -p ${DIR_OUT}
echo "BAM -> fastq for "${ID}
../scripts/bam2fastq.sh ${BAM_IN} ${BAM_IN}

echo "Now running bwa mem for "${ID}
../scripts/map2maizeAPGv4.sh ${ID} ${BAM_IN}_1.fq ${BAM_IN}_2.fq ${DIR_OUT}

echo "all done!"

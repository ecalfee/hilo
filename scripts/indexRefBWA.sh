#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J indexRefBWA
#SBATCH -o /home/ecalfee/hilo/slurm-log/indexRefBWA_%j.out
#SBATCH -t 10:00:00
#SBATCH -x bigmem1
#SBATCH --mem=30G

# general bash script settings to make sure if any errors in the pipeline fail
set –o pipefail
set –o errexit
set –o nounset

# load samtools
module load bwa
echo "indexing maize reference"
bwa index refMaize/AGPv4.fa
echo "all done!"

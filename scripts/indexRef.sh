#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /group/jrigrp/Share/assemblies/
#SBATCH -J indexRef
#SBATCH -o /home/ecalfee/hilo/slurm-log/indexRef_%j.out
#SBATCH -t 10:00:00
#SBATCH -x bigmem1
#SBATCH --mem=30G

# general bash script settings to make sure if any errors in the pipeline fail
set –o pipefail
set –o errexit
set –o nounset

# load samtools
module load samtools

samtools faidx Zea_mays.AGPv4.dna.chr.fa
echo "indexed maize reference"

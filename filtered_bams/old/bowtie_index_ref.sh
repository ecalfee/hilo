#!/bin/bash
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/data/refMaize
#SBATCH -J flagStat
#SBATCH -o /home/ecalfee/hilo/slurm-log/bowtieIndexRef_%A_%a.out
#SBATCH -t 6:00:00
#SBATCH --mem=24G

# to run: sbatch bowtie_index_ref.sh

# script creates bowtie2 index for reference genome (to use with fastq_screen)
# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load bio3
module load bowtie2

bowtie2-build -f Zea_mays.B73_RefGen_v4.dna.toplevel.fa Zea_mays.B73_RefGen_v4.dna.toplevel.fa

echo "all done!"

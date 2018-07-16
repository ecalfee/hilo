#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/
#SBATCH -J indexBam
#SBATCH -o /home/ecalfee/hilo/slurm-log/indexAllo4PopMaize_%j_%A_%a.out
#SBATCH -t 10:00:00
#SBATCH --mem=5G
#SBATCH --array=1-7,9-40

# general bash script settings to make sure if any errors in the pipeline fail
set –o pipefail
set –o errexit
set –o nounset

# set variables
# pad the task id with leading zeros
printf -v TASK_ID "%02g" $SLURM_ARRAY_TASK_ID
# necessary becuase slurm array task ID doesn't hold leading zeros

# load samtools
module load samtools

samtools index alloMaize4pop_symlink_bam/maizeLow_$TASK_ID.bam

echo "done indexing symlink "$TASK_ID".bam"

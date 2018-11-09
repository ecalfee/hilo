#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/scripts
#SBATCH -J geneDens
#SBATCH -o /home/ecalfee/hilo/slurm-log/getGeneDensity_%A_%a.out
#SBATCH -t 12:00:00
#SBATCH --mem=16G

# array contains all the chromosome #
CHR=$SLURM_ARRAY_TASK_ID

# directory for input and output is saved in the R file -- can change later

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load software needed
module load R


# run R script
echo "getting gene density for all sites"
Rscript rmap.R

echo 'all done!'

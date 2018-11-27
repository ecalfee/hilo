#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/scripts
#SBATCH -J xOgut
#SBATCH -o /home/ecalfee/hilo/slurm-log/extendOgut_%A_%a.out
#SBATCH -t 15:00
#SBATCH --mem=8G

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load software needed
module load R


# run R script
echo "extending Ogut recombination map to full range of chromosomes"
Rscript extendOgutMap2fullCHR.R
echo 'all done!'

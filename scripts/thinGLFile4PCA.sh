#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/scripts
#SBATCH -J thinGL
#SBATCH -o /home/ecalfee/hilo/slurm-log/thinGLFile4PCA_%A_%a.out
#SBATCH -t 15:00
#SBATCH --mem=8G
#SBATCH --array=0-425

# array contains all the region #'s'
REGION=$SLURM_ARRAY_TASK_ID

# directory for input and output
DIR_GL_IN="../data/geno_lik/merged_pass1_all_alloMaize4Low_16/allVar_depthFilt"
DIR_THINNED_SITES="../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedPCA"
DIR_OUT=${DIR_THINNED_SITES}"/GL_by_region"

mkdir -p ${DIR_OUT}

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load software needed
module load R


# run R script
echo "thinning GL for PCA: region "${REGION}
Rscript thinGLFile.R ${REGION} ${DIR_GL_IN} ${DIR_THINNED_SITES} ${DIR_OUT}

echo 'all done!'

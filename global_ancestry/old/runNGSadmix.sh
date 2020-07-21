#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/global_ancestry
#SBATCH -J NGSadmix
#SBATCH -o /home/ecalfee/hilo/slurm-log/NGSadmix_%A_%a.out
#SBATCH -t 6:00:00
#SBATCH --mem=8G
#SBATCH -n 1
#SBATCH --array=2-4

# to run: sbatch --export=PREFIX="hilo_alloMAIZE_MAIZE4LOW" runNGSadmix.sh

# slurm array task id sets number of genetic clusters, e.g.
# set an --array=2 for K = 2 or --array=2-4 to test K = 2, 3, 4 etc.

# set VARIABLES
k=$SLURM_ARRAY_TASK_ID
#DIR_GL="geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedPCA"
GL_FILE="results/thinnedSNPs/$PREFIX/whole_genome.beagle.gz"
OUT_DIR="results/NGSAdmix/$PREFIX"

# load module for NGSAdmix
module load bio

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# make directory to store output (if doesn't yet exist)
mkdir -p $OUT_DIR

echo "running NGSadmix"
NGSadmix -likes "${GL_FILE}" \
-K "${k}" -P 1 \
-o "${OUT_DIR}"/K"${k}"

# settings:
# -likes beagle genotype likelihood file
# -K 2 for number of subpopulations/clusters to consider in admixture model
# -P k splits the analysis job over k nodes (but does not distribute I/O)
# -o output

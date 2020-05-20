#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/global_ancestry
#SBATCH -J PCAngsd
#SBATCH -o /home/ecalfee/hilo/slurm-log/PCAngsd_%A_%a.out
#SBATCH -t 6:00:00
#SBATCH --mem=8G
#SBATCH -n 2

# to run: sbatch --export=PREFIX="hilo_alloMAIZE_MAIZE4LOW" runPCAngsd.sh

GL_FILE="results/thinnedSNPs/$PREFIX/whole_genome.beagle.gz"
DIR_OUT="results/PCA/$PREFIX"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load python 2 from bio module
module load bio #python/2.7.15

# run PCAngsd
mkdir -p "$DIR_OUT"

# assumes GL data is already filtered for a minimum MAF; doesn't re-filter
python2 /home/ecalfee/bin/pcangsd/pcangsd.py \
-beagle "$GL_FILE" \
-threads 2 -iter 100 \
-minMaf 0 -admix \
-o "$DIR_OUT/whole_genome"

# -admix option calculates admixture proportions in addition to genotype covariance matrix for PCA
# -iter specifies number of EM steps

echo "done running PCAngsd"

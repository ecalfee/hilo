#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J oultierPCAngsd
#SBATCH -o /home/ecalfee/hilo/slurm-log/outlierPCAngsd_%A_%a.out
#SBATCH -t 6:00:00
#SBATCH --mem=8G
#SBATCH -n 2

# note: must set the directory to print output and find regions.txt file
# when calling the script:
# to run:
# sbatch --export="DIR=outliers/chr4/inv4m,ALL" runPCAngsd_outliers.sh


GL_FILE_IN="$DIR/GL.beagle.gz"
PCA_FILE_OUT="$DIR/PCA"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load python 2 from bio module
module load bio

# run PCAngsd
echo "running PCAngsd for $GL_FILE_IN"

# assumes GL data is already filtered for a minimum MAF; doesn't re-filter
python2 /home/ecalfee/bin/pcangsd/pcangsd.py \
-beagle $GL_FILE_IN \
-threads 2 -iter 100 \
-minMaf 0 -admix \
-o $PCA_FILE_OUT

# -admix option calculates admixture proportions in addition to genotype covariance matrix for PCA
# -iter specifies number of EM steps

echo "done running PCAngsd"

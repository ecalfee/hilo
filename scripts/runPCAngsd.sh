#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J PCAngsd
#SBATCH -o /home/ecalfee/hilo/slurm-log/PCAngsd_%A_%a.out
#SBATCH -t 6:00:00
#SBATCH --mem=8G
#SBATCH -n 2

# the suffix for the genotype likelihood file in GL_PREFIX is omitted .beagle.gz;
GL_PREFIX="geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedPCA/whole_genome"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load python 2 from bio module
module load bio

# run PCAngsd
# assumes GL data is already filtered for a minimum MAF; doesn't re-filter
python2 /home/ecalfee/bin/pcangsd/pcangsd.py \
-beagle $GL_PREFIX.beagle.gz \
-threads 2 -iter 100 \
-minMaf 0 -admix \
-o $GL_PREFIX

# -admix option calculates admixture proportions in addition to genotype covariance matrix for PCA
# -iter specifies number of EM steps

echo "done running PCAngsd"

#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo
#SBATCH -J PCAngsd
#SBATCH -o /home/ecalfee/hilo/slurm-log/PCAngsd_%j_%A_%a.out
#SBATCH -t 36:00:00
#SBATCH --mem=30G
#SBATCH -n 10
#SBATCH --export=GL_PREFIX=data/geno_lik/pass1/allVar/whole_genome

# the suffix for the genotype likelihood file in GL_PREFIX is omitted .beagle.gz; 

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load python2 environment from anaconda
export PATH=/home/ecalfee/software/anaconda2/bin:$PATH
source activate condaEnv

# run PCAngsd
# assumes GL data is already filtered for a minimum MAF; doesn't re-filter
python2 ../bin/pcangsd/pcangsd.py -beagle $GL_PREFIX.beagle.gz -threads 10 -iter 100 -minMaf 0 -o $GL_PREFIX

echo "done running PCAngsd"

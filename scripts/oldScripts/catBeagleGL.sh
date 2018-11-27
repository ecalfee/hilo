#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/geno_lik/pass1/allVar
#SBATCH -J catGL
#SBATCH -o /home/ecalfee/hilo/slurm-log/catGL_%A_%a.out
#SBATCH -t 10:00:00
#SBATCH --mem=8G

# concatenates beagle genotype likelihood files, skipping headers

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# concatenate
zcat region_0.beagle.gz | gzip >> allVar.beagle.gz; for i in {1..46}; do zcat region_$i.beagle.gz | tail -n +2 | gzip >> whole_genome.beagle.gz; done

echo "done running with cat beagle.gz"

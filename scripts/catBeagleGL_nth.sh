#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/geno_lik/pass1/allVar
#SBATCH -J catGLNth
#SBATCH -o /home/ecalfee/hilo/slurm-log/catGLNth_%j_%A_%a.out
#SBATCH -t 10:00:00
#SBATCH --mem=8G

# concatenates beagle genotype likelihood files, skipping headers and only taking 1 in 100 positions

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# concatenate
zcat region_0.beagle.gz | awk 'NR % 100 == 0' | gzip >> whole_genome_pruned_by100.beagle.gz; for i in {1..46}; do zcat region_$i.beagle.gz | tail -n +2 | awk 'NR % 100 == 0' | gzip >> whole_genome_pruned_by100.beagle.gz; done

echo "done running with cat beagle.gz every 100"

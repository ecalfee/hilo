#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J NGSadmix
#SBATCH -o /home/ecalfee/hilo/slurm-log/NGSadmix_%j_%A_%a.out
#SBATCH -t 48:00:00
#SBATCH --mem=50G
#SBATCH -n 4

# slurm array task id sets number of genetic clusters, e.g.
# set an --array=2 for K = 2 or --array=2-4 to test K = 2, 3, 4 etc.

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# make directory to store output (if doesn't yet exist)
mkdir -p NGSadmix/pass1/

echo "running NGSadmix"
NGSadmix -likes geno_lik/pass1/pruned_all.beagle.gz \
-K $SLURM_ARRAY_TASK_ID -P 4 \
-o NGSadmix/pass1/K"$SLURM_ARRAY_TASK_ID"_pruned_all

# settings:
# -likes beagle genotype likelihood file
# -K 2 for number of subpopulations/clusters to consider in admixture model
# -P 4 splits the analysis job over 4 nodes (but does not distribute I/O)
# -o output
# NGSadmix is installed locally in user's bin & path

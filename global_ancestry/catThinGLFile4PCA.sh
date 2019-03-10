#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/global_ancestry
#SBATCH -J catThinGL
#SBATCH -o /home/ecalfee/hilo/slurm-log/catThinGL4PCA_%A.out
#SBATCH -t 4:00:00
#SBATCH --mem=8G

# concatenates 0-425 regions of beagle genotype likelihood files,
# skipping headers after first file (region 0)

# to run: sbatch --export=PREFIX=hilo_alloMAIZE_MAIZE4LOW catThinGL4PCA.sh


# set variables
DIR_OUT="results/thinnedSNPs/$PREFIX"
DIR_IN=${DIR_OUT}"/GL_by_region"
DIR_SCRATCH="/scratch/ecalfee/catThinGL4PCA/$PREFIX"
SCRATCH_OUT=${DIR_SCRATCH}"/whole_genome.beagle.gz"

# make local scratch DIRECTORY
mkdir -p "$DIR_SCRATCH"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# concatenate
zcat ${DIR_IN}/region_0.beagle.gz | gzip > ${SCRATCH_OUT}
for i in {1..425}
do zcat ${DIR_IN}/region_${i}.beagle.gz | tail -n +2 | gzip >> ${SCRATCH_OUT}
done

echo "done concatenating beagle.gz 4 PCA! Now transfering files from local to home directory"
rsync -avh --remove-source-files ${DIR_SCRATCH}/ ${DIR_OUT}/ # copies all contents of output directory over to the appropriate home directory & cleans up scratch dir.
echo "results copied to home output directory: "${DIR_OUT}

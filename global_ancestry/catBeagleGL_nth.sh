#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/global_ancestry
#SBATCH -J catGLNth
#SBATCH -o /home/ecalfee/hilo/slurm-log/catGLNth_%A_%a.out
#SBATCH -t 10:00:00
#SBATCH --mem=8G

# to run:
# sbatch --export=PREFIX=pass2pass2_alloMAIZE_PalmarChico,N=100 catBeagle_nth.sh

# concatenates beagle genotype likelihood files, skipping headers and only taking 1 in N positions
startR=0
endR=425
DIR_GL="../variant_sites/results/$PREFIX"
DIR_OUT="results/thinnedSNPs/$PREFIX/prunedBy$N"


# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# make output DIRECTORY
mkdir -p $DIR_OUT
# concatenate header from first file with every Nth SNP from each subsequent file
zcat $DIR_GL/region_$startR.beagle.gz | head -n 1 | gzip > $DIR_OUT/whole_genome.beagle.gz; \
for (( i=$startR; i<=$endR; i++ )); do zcat $DIR_GL/region_$i.beagle.gz | tail -n +2; done | awk -v N=$N 'NR % N == 0' | gzip >> $DIR_OUT/whole_genome.beagle.gz

echo "done running with cat beagle.gz every $N th position, regions $startN to $endN, sample set: $PREFIX"

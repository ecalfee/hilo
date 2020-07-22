#!/bin/bash

# to run:
# ./thin_GL_4PCA 100 0 425 HILO_MAIZE55

# concatenates beagle genotype likelihood files, skipping headers and only taking 1 in N positions
N=$1
startR=$2
endR=$3
PREFIX=$4
DIR_GL="variant_sites/results/$PREFIX"
DIR_OUT="global_ancestry/results/thinnedSNPs/$PREFIX"


# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# make output DIRECTORY
#mkdir -p $DIR_OUT
echo "working directory: "$PWD
echo "input: "$DIR_GL"/region_"$startR".beagle.gz"
echo "output:  "$DIR_OUT"/whole_genome.beagle.gz"
# concatenate header from first file with every Nth SNP from each subsequent file
zcat $DIR_GL/region_$startR.beagle.gz | head -n 1 | gzip > $DIR_OUT/whole_genome.beagle.gz; \
for (( i=$startR; i<=$endR; i++ )); do zcat $DIR_GL/region_$i.beagle.gz | tail -n +2; done | awk -v N=$N 'NR % N == 0' | gzip >> $DIR_OUT/whole_genome.beagle.gz

echo "done running with cat beagle.gz every $N th position, regions $startR to $endR, sample set: $PREFIX"

#!/bin/bash

# to run:
# ./concatenate_group_mafs.sh 0 425 "wavelets/results/alleleFreqs/HILO_MAIZE55/K2/byRegion/allopatric_maize wavelets/results/alleleFreqs/HILO_MAIZE55/K2/allopatric_maize.mafs.gz

# concatenates minor allele frequency files, skipping headers except for the first file
startR=$1
endR=$2
PREFIX_IN=$3
FILE_OUT=$4

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# make output DIRECTORY
#mkdir -p $DIR_OUT
echo "working directory: "$PWD
echo "input: $PREFIX_IN/region_$startR.mafs.gz"
echo "output: $FILE_OUT"
# concatenate header from first file with every Nth SNP from each subsequent file
zcat "$PREFIX_IN"/region_"$startR".mafs.gz | head -n 1 | gzip > "$FILE_OUT"; \
for (( i=$startR; i<=$endR; i++ )); do zcat "$PREFIX_IN"/region_$i.mafs.gz | tail -n +2; done | gzip >> "$FILE_OUT"

echo "all done!"

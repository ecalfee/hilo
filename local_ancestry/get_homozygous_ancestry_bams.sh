#!/bin/bash

# to run: local_ancestry/get_homozygous_ancestry_bams.sh {input.ploidy} {params.dir_tracts} {params.dir_in_bams} {params.dir_out}

# this script takes in a bed file for tracts with posterior probability > 0.8 of being homozygous for zea ancestry
# and filters all bams for a population to only include reads that overlap those regions


PLOIDY_FILE="$1"
DIR_TRACTS="$2"
DIR_IN_BAMS="$3"
DIR_OUT="$4"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset


echo $PLOIDY_FILE
echo $DIR_TRACTS
echo $DIR_IN_BAMS
echo $DIR_OUT

# make directory for output
mkdir -p "$DIR_OUT"

# get IDs from ploidy file
ids=($(cut -f1 $PLOIDY_FILE))
#ID=${ids["$i"]}

## how many ids in this pop?
n_ids=${#ids[@]}

#echo $n_ids


## Use bash for loop
for (( i=0; i<$n_ids; i++ )); do

		ID=${ids["$i"]}

		echo "filtering bam for homozygous ancestry regions $ID"

		bedtools intersect -sorted -a "$DIR_IN_BAMS/$ID.sort.dedup.bam" \
		-b "$DIR_TRACTS/$ID.bed" > "$DIR_OUT/$ID.sort.dedup.bam"

		echo "now indexing new bam!"
		sleep 5s # because index needs to have a later timestamp
		samtools index "$DIR_OUT/$ID.sort.dedup.bam"
done

echo "all done!"

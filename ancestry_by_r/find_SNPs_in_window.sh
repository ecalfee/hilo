#!/bin/bash

# gets genotype likelihoods from GL that overlap window in BED (sorted with scaffold order GENOME) -> sends to OUTPUT genotype likelihood file

# to run: ./find_SNPs_in_window.sh GL BED GENOME OUTPUT

GL=$1
BED=$2
GENOME=$3
OUTPUT=$4

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset


echo "GL: "$GL" BED: "$BED" GENOME: "$GENOME" OUTPUT: "$OUTPUT

zcat $GL | head -n 1 | gzip > $OUTPUT
zcat $GL | \
awk '{{split($1, a, "_"); print a[1]"\t"a[2]-1"\t"a[2]"\t"$0}}' | tail -n +2 | \
bedtools intersect -a stdin -b $BED -wa -sorted -g $GENOME | \
cut -f4- | gzip >> $OUTPUT

echo 'all done!'

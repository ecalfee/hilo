#!/bin/bash

# gets genotype likelihoods for block bootstrap sample (due to resampling, will have repeated SNPs)
# to run: ./make_bootstrap_GL_file.sh GL SAMPLE PREFIX OUTPUT

GL=$1
SAMPLE=$2
PREFIX=$3
OUTPUT=$4

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset


echo "GL: "$GL" SAMPLE: "$SAMPLE" PREFIX: "$PREFIX" OUTPUT: "$OUTPUT

(zcat $GL | head -n 1;
for w in $(cat $SAMPLE);
    do zcat ancestry_by_r/results/GL_1cM/$PREFIX/$w.beagle.gz | tail -n +2;
done) | gzip > $OUTPUT

echo 'all done!'

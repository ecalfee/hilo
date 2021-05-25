#!/bin/bash

# to run: get_homozygous_ancestry_tracts_K3.sh {input.bed} {params.ZEA} {input.ploidy} {params.dir_post} {params.dir_out}

# this script takes in a set of short tracts around sites with ancestry posteriors from ancestry_hmm
# and the MAP posteriors for individuals at those sites
# then filters for sites with posterior probability > 0.8 of being homozygous
# and merges tracts together with high probability homozygosity for zea ancestry

SITES_FILE="$1"
ZEA="$2"
PLOIDY_FILE="$3"
DIR_POST="$4"
DIR_OUT="$5"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset


echo $SITES_FILE
echo $ZEA
echo $PLOIDY_FILE
echo $DIR_POST
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

    if [ "$ZEA" == "maize" ]; then # find homozygous ancestry state for maize
        echo "finding maize tracts $ID"
		    # posterior probability for homozygous maize ancestry is the 8th column of ancestry_hmm posterior output file
        pr -mt -s$'\t' "$SITES_FILE" <(tail -n +2 "$DIR_POST/$ID.posterior" | cut -f8) | \
        awk '$5 > 0.8 {print $0}' | bedtools merge -sorted -d 1 > "$DIR_OUT/$ID.bed"

    elif [ "$ZEA" == "mexicana" ]; then
        echo "finding mexicana tracts $ID"
        # posterior probability for homozygous mexicana ancestry is the 3rd column of ancestry_hmm posterior output file
		    pr -mt -s$'\t' "$SITES_FILE" <(tail -n +2 "$DIR_POST/$ID.posterior" | cut -f3) | \
		    awk '$5 > 0.8 {print $0}' | bedtools merge -sorted -d 1 > "$DIR_OUT/$ID.bed"

		elif [ "$ZEA" == "parv" ]; then
			  echo "finding parviglumis (parv) tracts $ID"
	      # posterior probability for homozygous parviglumis ancestry is the 6th column of ancestry_hmm posterior output file
			  pr -mt -s$'\t' "$SITES_FILE" <(tail -n +2 "$DIR_POST/$ID.posterior" | cut -f6) | \
			  awk '$5 > 0.8 {print $0}' | bedtools merge -sorted -d 1 > "$DIR_OUT/$ID.bed"

    else
        echo "Input for ZEA needs to be maize or mexicana or parv!"

    fi

done

echo "all done!"

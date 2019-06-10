#!/bin/bash
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/within_ancestry_diversity
#SBATCH -J getMaize2
#SBATCH -o /home/ecalfee/hilo/slurm-log/getHomozygTractsBamsMaize_%A_%a.out
#SBATCH -t 2-00:00:00
#SBATCH --mem=16G
# to run: sbatch --array=0-217 --export=PREFIX=pass2_alloMAIZE,PREFIX_HMM=output_noBoot get_homozyg_tracts_and_bams.sh

# this script takes in a set of short tracts around sites with ancestry posteriors from ancestry_hmm
# and the posteriors for an individual at those sites
# then filters for sites with posterior probability > 0.8 of being homozygous
# and merges tracts together with high probability homozygosity for maize ancestry
# finally, filters bam files to only include aligned reads that overlap these homozygous ancestry tracts

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load bio # for bedtools and samtools

i=$SLURM_ARRAY_TASK_ID
ids=($(cat ../samples/${PREFIX}_IDs.list))
ID=${ids["$i"]}
bams=($(cat ../samples/${PREFIX}_bams.list))
BAM=${bams["$i"]}
SITES_FILE=results/input/"$PREFIX"/var.sites.bed
DIR_OUT="results/input/$PREFIX/mex2"
DIR_POST="../local_ancestry/results/ancestry_hmm/$PREFIX/$PREFIX_HMM"

# make directory for output
mkdir -p "$DIR_OUT/tracts"
mkdir -p "$DIR_OUT/bams"

# only run script if the posterior file can be found first
if [ ! -f "$DIR_POST"/"$ID".posterior ]; then
	echo "posterior does not exist: ${DIR_POST}/$ID.posterior"
else
  # getting high posterior probability tracts homozygous for mexicana ancestry
  echo "finding mexicana tracts"
  pr -mt -s$'\t' "$SITES_FILE" <(tail -n +2 "$DIR_POST/$ID".posterior | cut -f3) | \
  awk '$4 > 0.8 {print $0}' | bedtools merge -sorted -d 1 > $DIR_OUT/tracts/"$ID".bed

  # filtering bam for homozygous mex regions
  echo "filtering bam for homozygous mexicana regions"
  bedtools intersect -sorted -a "$BAM" \
  -b $DIR_OUT/tracts/"$ID".bed > $DIR_OUT/bams/"$ID".sort.dedup.baq.bam

  echo "now indexing new bam!"
  sleep 5s # because index needs to have a later timestamp
  samtools index $DIR_OUT/bams/"$ID".sort.dedup.baq.bam

  echo "all done!"
fi

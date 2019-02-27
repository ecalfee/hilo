#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/within_ancestry_diversity
#SBATCH -J getMex2
#SBATCH -o /home/ecalfee/hilo/slurm-log/getHomozygTractsBamsMex_%A_%a.out
#SBATCH -t 2-00:00:00
#SBATCH --mem=16G
#SBATCH --array=1-200
# to run: sbatch --export=DIR_BAMS=../filtered_bams/results/pass1,DIR_POST=../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/output_noBoot get_homozyg_tracts_and_bams.sh

# this script takes in a set of short tracts around sites with ancestry posteriors from ancestry_hmm
# and the posteriors for an individual at those sites
# then filters for sites with posterior probability > 0.8 of being homozygous
# and merges tracts together with high probability homozygosity formexicana ancestry
# finally, filters bam files to only include aligned reads that overlap these homozygous ancestry tracts

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load bio # for bedtools and samtools

i=$SLURM_ARRAY_TASK_ID
SITES_FILE=results/input/var.sites.bed

# make directory for output
mkdir -p "results/input/mex2/tracts"
mkdir -p "results/input/mex2/bams"

# only run script if the posterior file can be found first
if [ ! -f "$DIR_POST"/HILO"$i".posterior ]; then
	echo "posterior does not exist: ${DIR_POST}/HILO${i}.posterior"
else
  # getting high posterior probability tracts homozygous for mexicana ancestry
  echo "finding mexicana tracts"
  pr -mt -s$'\t' "$SITES_FILE" <(tail -n +2 "$DIR_POST"/HILO"$i".posterior | cut -f5) | \
  awk '$4 > 0.8 {print $0}' | bedtools merge -sorted -d 1 > results/input/mex2/tracts/HILO"$i".bed

  # filtering bam for homozygous mex regions
  echo "filtering bam for homozygous mexicana regions"
  bedtools intersect -sorted -a "$DIR_BAMS"/HILO"$i".sort.dedup.baq.bam \
  -b results/input/mex2/tracts/HILO"$i".bed > results/input/mex2/bams/HILO"$i".sort.dedup.baq.bam

  echo "now indexing new bam!"
  sleep 5s # because index needs to have a later timestamp
  samtools index results/input/mex2/bams/HILO"$i".sort.dedup.baq.bam

  echo "all done!"
fi

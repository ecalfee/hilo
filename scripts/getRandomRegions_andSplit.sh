#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J splitRand
#SBATCH -o /home/ecalfee/hilo/slurm-log/getRandomRegionsandSplit_%j_%A_%a.out
#SBATCH -t 5:00:00
#SBATCH --mem=2G
#SBATCH --export=N=22000,L=100,SEED=712,LINES_PER_FILE=500

# default is to generate 22,000 regions length 100 bp with random seed 712 and split them into (44 files) with 500 regions each
# this covers ~ 1/1000 of the maize genome

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load software needed
module load bedtools

# make directory to store output (if doesn't yet exist)
mkdir -p refMaize/random_regions/N$N.L$L.regions_split

# generate regions and sort by chrom and position within chrom and then format chr:start-end
bedtools random -g refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths -n $N -l $L -seed $SEED\
| sort -nbk1,1 -nbk2,2 | awk '{print $1"\:"$2"-"$3}' > refMaize/random_regions/N$N.L$L.regions

# split; d is for labelling split files with digits
split -d -l $LINES_PER_FILE refMaize/random_regions/N$N.L$L.regions refMaize/random_regions/N$N.L$L.regions_split/split_

echo 'done printing and splitting regions'

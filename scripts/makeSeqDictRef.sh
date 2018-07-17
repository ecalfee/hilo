#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/refMaize
#SBATCH -J dictRef
#SBATCH -o /home/ecalfee/hilo/slurm-log/makeSeqDictRef_%j.out
#SBATCH -t 16:00:00
#SBATCH --mem=12G

# this script makes symlinks to the shared reference genome AGPv4 and saves in my local directory
# then it makes a sequence dictionary for the reference, as required by picard

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# make symlinks to ref.fa and index ref.fa.fai
ln -s /group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa AGPv4.fa
ln -s /group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa.fai AGPv4.fa.fai

# load java to run picardtools
module load java
# load picardtools for removing PCR duplicate reads
module load picardtools # saves path to loaded versin in $PICARD variable

# create sequence dictionary for AGP4 reference
java -jar $PICARD/picard.jar CreateSequenceDictionary \
      R=AGPv4.fa \
      O=AGPv4.dict

echo "done creating symlinks and sequence dictionary for APGv4"

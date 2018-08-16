#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J SAFAngsd
#SBATCH -o /home/ecalfee/hilo/slurm-log/SAFAngsd_%j_%A_%a.out
#SBATCH -t 3:00:00
#SBATCH -n 1
#SBATCH --mem=2G
#SBATCH --export=REGION_FILE=N1000.L100.regions
#SBATCH --array=18-31,33-35,360-363,365-374,1000,2000,3000

# pops >=1000 are the larger groupings of maize sympatric, mexicana sympatric and mexicana allopatric inds

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# rather than load module, use local module of angsd (v9.20)

# make directory to store output (if doesn't yet exist)
mkdir -p SAF/pass1/$REGION_FILE

echo "finding allele frequencies for pop"$SLURM_ARRAY_TASK_ID", regions"$REGION_FILE
# steps:
# (0) Start with filtered BAM files and reference genome
# (1) Use all sites to estimate site allele frequency
angsd -out SAF/pass1/$REGION_FILE/pop$SLURM_ARRAY_TASK_ID \
-anc refMaize/AGPv4.fa -fold 1 \
-rf refMaize/random_regions/$REGION_FILE \
-bam pass1_bam_pops/pop$SLURM_ARRAY_TASK_ID.list \
-remove_bads 1 -minMapQ 30 \
-GL 1 -dosaf 1 \
-P 1
echo "done with SAF pop"$SLURM_ARRAY_TASK_ID" regions "$REGION_FILE

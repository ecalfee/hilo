#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J SFS1DAngsd
#SBATCH -o /home/ecalfee/hilo/slurm-log/SFS1DAngsd_%j_%A_%a.out
#SBATCH -t 2:00:00
#SBATCH --mem=2G
#SBATCH --export=REGION_FILE=N1000.L100.regions
#SBATCH --array=18-31,33-35,360-363,365-374,1000,2000,3000,5000

# pops >=1000 are the larger groupings of maize sympatric, mexicana sympatric and mexicana allopatric inds

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

POP=$SLURM_ARRAY_TASK_ID

# rather than load module, use local module of angsd (v9.20)

# make directory to store output (if doesn't yet exist)
mkdir -p SFS/pass1/$REGION_FILE

echo "finding SFS for pop"$SLURM_ARRAY_TASK_ID", regions"$REGION_FILE
# steps:
# (0) Start with filtered BAM files and reference genome
# (1) Use all sites to estimate the folded site allele frequency using script calcSAFAngsd.sh
# (2) Then get the maximum likelihood estimate of the SFS
realSFS SAF/pass1/$REGION_FILE/pop$POP.saf.idx -P 1 > SFS/pass1/$REGION_FILE/pop$POP.folded.sfs

# (3) Calculate thetas -- actually not needed right now for whole-genome approach (vs. scanning for selection/thetas)
# I will instead calculate pi directly from the estimated genome-wide SFS
# for scanning approach: http://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests
#angsd -bam pass1_bam_pops/pop$POP.list \
#-out SFS/pass1/$REGION_FILE/pop$POP \
#-rf refMaize/random_regions/$REGION_FILE -doThetas 1 \
#-doSaf 1 -pest SFS/pass1/$REGION_FILE/pop$POP.folded.sfs \
#-anc refMaize/AGPv4.fa -GL 1 -fold 1 \
#-remove_bads 1 -minMapQ 30
# -rf indicates to analyze a subset of regions from the whole genome
# uses maximum likelihood SFS as a prior 


#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J calcDepthRegions
#SBATCH -o /home/ecalfee/hilo/slurm-log/calcDepthRegions_%A_%a.out
#SBATCH -t 36:00:00
#SBATCH --mem=10G
#SBATCH --export=REGIONS_FILE=N1000.L100.regions

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# make directory to store output (if doesn't yet exist)
mkdir -p depthCov/pass1

# apply filtering with SAMtools & PICARD
echo "calculating depth using ANGSD on BAMS for hilo "$REGIONS_FILE
# calculate depth only including bases meeting minQ >20
angsd -out depthCov/pass1/$REGIONS_FILE.Q20 \
-rf refMaize/random_regions/$REGIONS_FILE \
-bam pass1_bam.all.list \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doCounts 1 \
-doDepth 1 
# calculate depth for bases meeting any base quality
angsd -out depthCov/pass1/$REGIONS_FILE \
-rf refMaize/random_regions/$REGIONS_FILE \
-bam pass1_bam.all.list \
-remove_bads 1 \
-minMapQ 30 \
-doCounts 1 \
-doDepth 1 
echo "done calculating depth using ANGSD on BAMS for hilo "$REGIONS_FILE

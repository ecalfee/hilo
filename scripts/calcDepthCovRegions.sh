#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J calcDepthRegions
#SBATCH -o /home/ecalfee/hilo/slurm-log/calcDepthRegions_%A_%a.out
#SBATCH -t 4:00:00
#SBATCH --mem=8G

# set some variables
REGIONS_LIST=N1000.L100.regions
OUT_DIR=depthCov/merged_pass1_all_alloMaize4Low_16
BAM_LIST=merged_bam.pass1_all.alloMaize4Low_16.list

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# make directory to store output (if doesn't yet exist)
mkdir -p ${OUT_DIR}

# apply filtering with SAMtools & PICARD
echo "calculating depth using ANGSD on BAMS for hilo "$REGIONS_FILE
# calculate depth only including bases meeting minQ >20
angsd -out ${OUT_DIR}/${REGIONS_LIST}.Q20 \
-rf refMaize/random_regions/${REGIONS_LIST} \
-bam ${BAM_LIST} \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doCounts 1 \
-doDepth 1 
echo "done calculating depth base quality > 20 & now calculating across all reads that map"
# calculate depth for bases meeting any base quality
angsd -out ${OUT_DIR}/${REGIONS_LIST} \
-rf refMaize/random_regions/${REGIONS_LIST} \
-bam ${BAM_LIST} \
-remove_bads 1 \
-minMapQ 30 \
-doCounts 1 \
-doDepth 1 
echo "done calculating depth using ANGSD on BAMS for hilo "$REGIONS

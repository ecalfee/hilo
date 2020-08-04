#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/variant_sites
#SBATCH -J calcDepthRegions
#SBATCH -o /home/ecalfee/hilo/slurm-log/calcDepthRegions_%A.out
#SBATCH -t 4:00:00
#SBATCH --mem=8G

# to run: sbatch --export=PREFIX=hilo,REGIONS_LIST=N1000.L100.regions calcDepthCovRegions.sh

# set some variables
OUT_DIR=results/depthCov/"$REGIONS_LIST"
BAM_LIST=../samples/"$PREFIX"_bams.list
MAX_DEPTH=10000
# MAX_DEPTH is the maximum number of reads counted per site -- any more is put in a bin >=MAX_DEPTH

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load bio #angsd

# make directory to store output (if doesn't yet exist)
mkdir -p ${OUT_DIR}

# apply filtering with SAMtools & PICARD
echo "calculating depth using ANGSD on BAMS for hilo "$REGIONS_FILE
# calculate depth only including bases meeting minQ >20
angsd -out ${OUT_DIR}/${PREFIX}.Q20 \
-rf ../data/refMaize/random_regions/${REGIONS_LIST} \
-bam ${BAM_LIST} \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doCounts 1 \
-doDepth 1 \
-maxDepth ${MAX_DEPTH}
echo "done calculating depth base quality > 20 & now calculating across all reads that map"
# calculate depth for bases meeting any base quality
angsd -out ${OUT_DIR}/${PREFIX} \
-rf ../data/refMaize/random_regions/${REGIONS_LIST} \
-bam ${BAM_LIST} \
-remove_bads 1 \
-minMapQ 30 \
-doCounts 1 \
-doDepth 1 \
-maxDepth ${MAX_DEPTH}
echo "done calculating depth using ANGSD on BAMS for hilo "$REGIONS

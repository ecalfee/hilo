#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/local_ancestry
#SBATCH -J prune4HMM
#SBATCH -o /home/ecalfee/hilo/slurm-log/pruneFixcM4HMM_%A_%a.out
#SBATCH -t 10:00:00
#SBATCH --mem=8G
#SBATCH --array=1-10

# to run: sbatch --export=PREFIX=hilo_alloMAIZE_MAIZE4LOW pruneFixedcM_HMM.sh

# this script prunes variant sites to ~unlinked fixed distance (.001cM~1kb)
# for analysis with global ancestry and PCA
i=$SLURM_ARRAY_TASK_ID #chromosome
REGIONS_LIST=../data/refMaize/divide_5Mb/ALL_regions.list
MIN_cM=0.001
DIR_SITES="results/ancestryInformativeSNPs/$PREFIX"
DIR_OUT="results/thinnedSNPs/$PREFIX"
DIR_SCRATCH="/scratch/ecalfee/thinSNPs4HMM_chr"${i}

# general bash script settings to make sure if any errors in the pipeline fail
set –o pipefail
set –o errexit
set –o nounset

# load modules
# using python3 for pandas
# and bio for angsd v9.21 for angsd sites indexing
#module load python3
module load bio3
# load module bio later in script to not conflict with python3 in bio3

# make output directory
mkdir -p ${DIR_OUT}
mkdir -p ${DIR_SCRATCH}

echo "pruning SNPs"
python3 ../variant_sites/pruneFixedcM.py ${MIN_cM} ${DIR_SCRATCH}"/chr"${i} \
$(for i in $(awk -v chr=${i} '$1 == chr {print $4}' ${REGIONS_LIST}); \
do echo ${DIR_SITES}/region_${i}.var.sites; done)
# second line finds all regions associated with a specific chromosome (from file ALL_regions.list) and lists those mafs.gz files as input files to pruneFixedcM.py
echo "done pruning SNPs chr"$i

echo "fixing end of line character"
# copy to temporary file name then delete extra newline character when copied back to original file name
cp ${DIR_SCRATCH}/chr${i}.var.sites ${DIR_SCRATCH}/temp_chr${i}.var.sites
cat ${DIR_SCRATCH}/temp_chr${i}.var.sites | tr -d '\r' > ${DIR_SCRATCH}/chr${i}.var.sites
# remove temporary sites file
rm ${DIR_SCRATCH}/temp_chr${i}.var.sites

module load bio # for angsd

echo "making ANGSD sites index"
sleep 2s # wait 2 seconds before indexing so that index doesn't have same timestamp as sites file
angsd sites index ${DIR_SCRATCH}/chr${i}.var.sites
echo "all done!"

echo "all done making thinned sites! Now transfering files from local to home directory"
rsync -avh --remove-source-files ${DIR_SCRATCH}/ ${DIR_OUT}/ # copies all contents of output directory over to the appropriate home directory & cleans up scratch dir.
echo "results copied to home output directory: "${DIR_OUT}

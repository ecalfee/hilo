#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/
#SBATCH -J pruneFixcM
#SBATCH -o /home/ecalfee/hilo/slurm-log/pruneFixcM_%A_%a.out
#SBATCH -t 10:00:00
#SBATCH --mem=8G
#SBATCH --array=1-10
#SBATCH --export=MIN_cM=0.001,dir_in="var_sites/merged_pass1_all_alloMaize4Low_16/filteredAlloMAFnInd",dir_out="var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM"

# this script prunes variant sites to fixed distance (.001cM~1kb)
i=$SLURM_ARRAY_TASK_ID #chromosome
regions_list=refMaize/divide_5Mb/ALL_regions.list

# general bash script settings to make sure if any errors in the pipeline fail
set –o pipefail
set –o errexit
set –o nounset

# load modules
# using python3 for pandas
module load python3
# do not load ANGSD (local version is most up-to-date v9.20)

# make output directory
mkdir -p ${dir_out}

echo "pruning SNPs"
python3 ../scripts/pruneFixedcM.py ${MIN_cM} ${dir_out}"/chr"${i} \
$(for i in $(awk -v chr=${i} '$1 == chr {print $4}' ${regions_list}); \
do echo ${dir_in}/region_${i}.var.sites; done)
# second line finds all regions associated with a specific chromosome (from file ALL_regions.list) and lists those mafs.gz files as input files to pruneFixedcM.py
echo "done pruning SNPs chr"$i

echo "fixing end of line character"
# copy to temporary file name then delete extra newline character when copied back to original file name
cp ${dir_out}/chr${i}.var.sites ${dir_out}/temp_chr${i}.var.sites
cat ${dir_out}/temp_chr${i}.var.sites | tr -d '\r' > ${dir_out}/chr${i}.var.sites
# remove temporary sites file
rm ${dir_out}/temp_chr${i}.var.sites

echo "making ANGSD sites index"
sleep 2s # wait 2 seconds before indexing so that index doesn't have same timestamp as sites file
angsd sites index ${dir_out}/chr${i}.var.sites
echo "all done!"

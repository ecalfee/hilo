#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J countACGT
#SBATCH -o /home/ecalfee/hilo/slurm-log/countACGT_%A_%a.out
#SBATCH -t 4:00:00
#SBATCH --mem=8G
#SBATCH --array=1-200
#SBATCH --export="CHR=1,ALL"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# this script takes in a sites file and outputs an angsd counts file with
# columns totA totC totG totT as well as a .pos.gz file with all positions with data
# for a low-coverage hilo individual

# array is the hilo ID #
id=$SLURM_ARRAY_TASK_ID
dir="geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM"
bam_dir="filtered_bam"
bam_path=${bam_dir}"/hilo_"${id}".sort.dedup.baq.bam"
DIR_OUT=${dir}"/countsACGT/chr"${CHR}
DIR_SCRATCH="/scratch/ecalfee/countACGT_hilo_"${id}"_chr"${i}
MAX_DEPTH_IND=8

# load modules
# includes angsd v9.21
module load bio

# make output and scratch directories
mkdir -p ${DIR_OUT}
mkdir -p ${DIR_SCRATCH}

# make file with single focus bam listed
mkdir -p ${bam_dir}/list_indiv_bams
echo ${bam_path} > ${bam_dir}/list_indiv_bams/hilo_${id}.only.list

echo "counting reads supporting ACGT alleles for chr"${CHR}" and individual HILO"${id}
# counts reads for sites in sites file that pass basic quality filtering
angsd -out ${dir}/countsACGT/chr${CHR}/hilo_${id} \
-ref refMaize/AGPv4.fa \
-bam ${bam_dir}/list_indiv_bams/hilo_${id}.only.list \
-minQ 20 -minMapQ 30 \
-doCounts 1 -dumpCounts 3 \
-setMaxDepthInd ${MAX_DEPTH_IND} \
-remove_bads 1 \
-r ${CHR} \
-sites ${dir}/chr${CHR}.var.sites

# transfer output files to home directory
echo "all done counting ACGT! Now transfering files from local to home directory"
rsync -avh --remove-source-files ${DIR_SCRATCH}/ ${DIR_OUT}/ # copies all contents of output directory over to the appropriate home directory & cleans up scratch dir.
echo "results copied to home output directory: "${DIR_OUT}

echo "all done!"

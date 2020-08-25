#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/local_ancestry
#SBATCH -J ancInform
#SBATCH -o /home/ecalfee/hilo/slurm-log/filterAncestryInformativeAlloMAFnInd_%A_%a.out
#SBATCH -t 1:00:00
#SBATCH --mem=8G
#SBATCH --array=0-425

# to run: sbatch --export=PREFIX=pass2_alloMAIZE,D=0.3,MIN_MAIZE=10,MIN_MEX=10 filterAlloMAFnInd.sh

# other variables
#D=0.3
#MIN_MAIZE=10
#MIN_MEX=10
allo_maize="allo.maize"
allo_mex="allo.mexicana"
dir_maf="../within_ancestry_diversity/results/allele_freq/"$PREFIX
dir_sites="../variant_sites/results/"$PREFIX
dir_out="results/ancestryInformativeSNPs/"$PREFIX


region=$SLURM_ARRAY_TASK_ID

# general bash script settings to make sure if any errors in the pipeline fail
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load R

# make output directory
mkdir -p $dir_out

# filter for informative sites by minimum allele freq difference and min. # individuals with data:
# Rscript args (in order): region, MIN_MAIZE, MIN_MEX, D, allo_maize, allo_mex, dir_gl, dir_sites, dir_out
Rscript ./filter_var_sites_min_allo.R $region $MIN_MAIZE $MIN_MEX $D $allo_maize $allo_mex \
$dir_maf $dir_sites $dir_out

echo "done filtering for informative SNPs in R for region "$region

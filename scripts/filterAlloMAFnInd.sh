#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/
#SBATCH -J filterAllo
#SBATCH -o /home/ecalfee/hilo/slurm-log/filterAlloMAFnInd_%j_%A_%a.out
#SBATCH -t 1:00:00
#SBATCH --mem=8G
#SBATCH --array=0-425
#SBATCH --export=minInd=6,D=0.3

# other variables
allo_maize="maize.allo.4Low16"
allo_mex="mexicana.allo.withXochi35"
dir_gl="geno_lik/merged_pass1_all_alloMaize4Low_16/allVar/"
dir_sites="var_sites/merged_pass1_all_alloMaize4Low_16/"
dir_out="var_sites/merged_pass1_all_alloMaize4Low_16/filteredAlloMAFnInd"

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
# Rscript args (in order): region, minInd, D, allo_maize, allo_mex, dir_gl, dir_sites, dir_out
Rscript ../scripts/filter_var_sites_min_allo.R $region $minInd $D $allo_maize $allo_mex \
$dir_gl $dir_sites $dir_out

echo "done filtering for informative SNPs in R for region "$region


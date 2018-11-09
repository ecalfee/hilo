#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/
#SBATCH -J ancInform
#SBATCH -o /home/ecalfee/hilo/slurm-log/filterAncestryInformativeAlloMAFnInd_%A_%a.out
#SBATCH -t 1:00:00
#SBATCH --mem=8G
#SBATCH --array=0-425

# other variables
D=0.3
minInd_maize=10
minInd_mex=4
allo_maize="maize.allo.4Low16"
allo_mex="mexicana.allo.withXochi35"
dir_main="geno_lik/merged_pass1_all_alloMaize4Low_16"
dir_maf=${dir_main}/"allVar_depthFilt"
dir_sites=${dir_maf}
dir_out=${dir_main}"/ancestryInformativeFilt"

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
# Rscript args (in order): region, minInd_maize, minInd_mex, D, allo_maize, allo_mex, dir_gl, dir_sites, dir_out
Rscript ../scripts/filter_var_sites_min_allo.R $region $minInd_maize $minInd_mex $D $allo_maize $allo_mex \
$dir_maf $dir_sites $dir_out

echo "done filtering for informative SNPs in R for region "$region

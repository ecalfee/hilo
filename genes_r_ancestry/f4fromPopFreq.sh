#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/genes_r_ancestry
#SBATCH -J f4
#SBATCH -o /home/ecalfee/hilo/slurm-log/f4fromPopFreq_%A.out
#SBATCH -t 6:00:00
#SBATCH --mem=8G
#SBATCH --export="POP1=parv,POP2=symp.maize,POP3=allo.mexicana,POP4=trip,PREFIX=all"


# directory for input and output
DIR_IN="../within_ancestry_diversity/results/allele_freq/"$PREFIX
DIR_OUT="results/f4/"$PREFIX"/"$POP1"_"$POP2"_"$POP3"_"$POP4

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load software needed
module load R

# make output directory as needed
mkdir -p "$DIR_OUT"

# run R script
echo "calculating f4 statistic"
Rscript ./f4_from_pop_freq.R ${POP1} ${POP2} ${POP3} ${POP4} ${DIR_IN} ${DIR_OUT}

echo 'all done!'

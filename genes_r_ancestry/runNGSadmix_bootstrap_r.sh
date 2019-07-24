#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/genes_r_ancestry
#SBATCH -J NGSadmix_boot
#SBATCH -o /home/ecalfee/hilo/slurm-log/NGSadmix_bootstrap_r_%A_%a.out
#SBATCH -t 12:00:00
#SBATCH --mem=8G
#SBATCH -n 1
#SBATCH --array=0-100

# to run: sbatch --export=PREFIX=pass2_alloMAIZE_PalmarChico,K=3,WIND=1cM runNGSadmix_bootstrap_r.sh

# slurm array task id sets bootstrap.
# K sets number of genetic clusters, e.g.
# set an --array=2 for K = 2 or --array=2-4 to test K = 2, 3, 4 etc.
# r sets the recombination bin (out of 5)

# set VARIABLES
BOOT=$SLURM_ARRAY_TASK_ID

# load module for NGSAdmix
module load bio

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

echo "running NGSadmix"
for r in {1..5}
  do echo "recombination bin $r"
  DIR="results/bootstrap/windows_$WIND/r5_recomb$r/$PREFIX"
  GL_FILE="$DIR/boot$BOOT.beagle.gz"
  mkdir -p "$DIR/K${K}"
  NGSadmix -likes "${GL_FILE}" \
  -K "${K}" -P 1 \
  -o "${DIR}"/K"${K}"/boot"${BOOT}"
  echo "created output file: ${DIR}/K${K}/boot${BOOT}"
  echo "using R to define which ancestry is mexicana"
  Rscript ./define_K_ancestries.R ${PREFIX} ${r} ${WIND} ${BOOT} ${K}
done
echo "all done!"
# settings:
# -likes beagle genotype likelihood file
# -K 2 for number of subpopulations/clusters to consider in admixture model
# -P k splits the analysis job over k nodes (but does not distribute I/O)
# -o output

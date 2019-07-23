#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/genes_r_ancestry
#SBATCH -J bootstrapGL
#SBATCH -o /home/ecalfee/hilo/slurm-log/makeBootstrapGLFiles_%A_%a.out
#SBATCH -t 1:00
#SBATCH --mem=8G
#SBATCH --array=1-100

# to run: genes_r_ancestry$ sbatch --array=1-100 --export=PREFIX=pass2_alloMAIZE_PalmarChico,r_bins=5,prunedBy=100,WIND=1cM makeBootstrapGLFiles.sh

i=$SLURM_ARRAY_TASK_ID # index of bootstrap

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# run shell command
for r in {1..r_bins};
  do (zcat ../global_ancestry/results/thinnedSNPs/$PREFIX/prunedBy$prunedBy/whole_genome.beagle.gz | head -n 1
  for w in $(cat results/bootstrap/windows_$WIND/r5_recomb$r/boot$i.list)
    # make output DIRECTORY
    DIR_OUT=results/bootstrap/windows_$WIND/r${r_bins}_recomb$r/$PREFIX
    mkdir -p $DIR_OUT
    do zcat results/GL_$WIND/$PREFIX/$w.beagle.gz; done) | \
    gzip > $DIR_OUT/boot$i.beagle.gz;
  done;
done
echo 'all done!'

#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J catGLNth
#SBATCH -o /home/ecalfee/hilo/slurm-log/catGLNth_%j_%A_%a.out
#SBATCH -t 10:00:00
#SBATCH --mem=8G
#SBATCH --export=startR=0,endR=46,DIR_GL=data/geno_lik/pass1/allVar,N=1000

# concatenates beagle genotype likelihood files, skipping headers and only taking 1 in N positions

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# concatenate header from first file with every Nth SNP from each subsequent file
zcat $DIR_GL/region_$startR.beagle.gz | head -n 1 | gzip > $DIR_GL/whole_genome_pruned_every_$N.beagle.gz; \
for i in {$startR..$endR}; do zcat region_$i.beagle.gz | tail -n +2; done | awk -v N=$N 'NR % N == 0' | gzip >> $DIR_GL/whole_genome_pruned_every_$N.beagle.gz

echo "done running with cat beagle.gz every "$N"th position, regions "$startN" to "$endN

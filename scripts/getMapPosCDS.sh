#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/
#SBATCH -J mapCDS
#SBATCH -o /home/ecalfee/hilo/slurm-log/mapCDS_%A_%a.out
#SBATCH -t 4:00:00
#SBATCH --mem=8G
#SBATCH --array=1-10

# this script uses a python script to get map positions for CoDing Sequence ranges
i=$SLURM_ARRAY_TASK_ID #chromosome
fileIn="refMaize/geneAnnotations/CDS_"$i".txt"
fileOut="refMaize/geneAnnotations/CDS_"$i"_mapPos.txt"

# general bash script settings to make sure if any errors in the pipeline fail
set –o pipefail
set –o errexit
set –o nounset

# load modules
# using python3 for pandas
module load python3

echo "getting map positions"
python3 ../scripts/getMapPosCDS.py ${i} ${fileIn} ${fileOut}

echo "all done!"

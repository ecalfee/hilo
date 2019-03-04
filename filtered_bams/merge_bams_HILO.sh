#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/filtered_bams
#SBATCH -J mergeBam
#SBATCH -o /home/ecalfee/hilo/slurm-log/mergeBamsHILO_%A_%a.out
#SBATCH -t 1-00:00:00
#SBATCH --mem=32G
#SBATCH -n 1
#SBATCH --array=1-267

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
# load bio for samtools and bwa
module load bio
# load java to run picardtools
module load java
# load picardtools for removing PCR duplicate reads
module load picardtools # saves path to loaded versin in $PICARD variable

echo "working directory: ${PWD}" # print current working directory
echo "picard path: ${PICARD}"

ID=HILO"$SLURM_ARRAY_TASK_ID"
DIR_OUT=results/merged
DIR_METRICS=metrics/merged # metrics directory
DIR_TMP="/scratch/ecalfee/${ID}" # for memory overflow

# make output directory
mkdir -p ${DIR_OUT}
mkdir -p ${DIR_METRICS}
mkdir -p ${DIR_TMP}


# make results directory (if necessary)
mkdir -p merged
echo "working directory:"${PWD} # print current working directory

echo "merging sorted bams"

# are there different files to merge?
if [ $(wc -l results/to_merge/"$ID".list | cut -f1 -d" ") = 1 ]; then
	echo "no files to merge for $ID. creating symlink"
	ln -s $(cat results/to_merge/"$ID".list) "$DIR_OUT"/"$ID".sort.dedup.baq.bam;
	ln -s $(cat results/to_merge/"$ID".list)".bai" "$DIR_OUT"/"$ID".sort.dedup.baq.bam.bai;
else
	echo "merging files:"
	cat results/to_merge/"$ID".list
	
	# merging
	samtools merge -b results/to_merge/"$ID".list "$DIR_OUT"/"$ID".sort.baq.bam

	# need to deduplicate again because even though sequenced on diff. lanes they are the same libraries,
	# so share PCR duplicates. I don't need to redo sorting or BAQ scores however.
	
	echo "marking duplicates with PICARD"
	java -Xmx6g -jar ${PICARD}/picard.jar MarkDuplicates \
	INPUT="${DIR_OUT}/${ID}.sort.baq.bam" OUTPUT=/dev/stdout QUIET=true \
	REMOVE_DUPLICATES=true \
	TMP_DIR="${DIR_TMP}" \
	METRICS_FILE="${DIR_METRICS}/${ID}.metrics.txt" |\
	samtools view -bS - > "${DIR_OUT}/${ID}.sort.dedup.baq.bam"

	echo "removing intermediate bam"
	rm "${DIR_OUT}/${ID}.sort.baq.bam"

	echo "now indexing merged and deduplicated bam!"
	sleep 5s # because index needs to have a later timestamp
	samtools index "${DIR_OUT}/${ID}.sort.dedup.baq.bam"

	echo "all done!"
fi

# options

# picard markduplicates:
# -Xmx6g will not spill to tmp until 6G of memory are used

# samtools converts sam to bam output:
# -S specifies that input is sam, not binary bam
# - sets stdin as input (stdout as output is default)
# -b is for bam output

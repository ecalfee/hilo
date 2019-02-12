After read mapping, we use the whole sample (reference and admixed individuals) to identify variant sites using ANGSD. Steps:
1) Identify sites with >5% MAF in global sample based on reads meeting mapping and base quality cutoffs. 
Filters out sites with too high coverage in global sample and sites with data for too few individuals.
mafs.gz output files from angsd are divided into 5Mb regions across the genome, regions 0-425:
allVarAngsdGL1.sh
(note: high coverage cutoff was based on calculating total depth of coverage (from ANGSD) across some random regions: 
R script plot_mapping_metrics.R was used to process the results of calcDepthCovRegions.sh)

2) Create and index 'sites' file for each mafs.gz file, listing the chromosome and bp position of each variant site:
mafsToSitesFile.sh

Note: Filtering of individual level data for low or high coverage happens within each analysis.
Different analyses require different filtering from this point for levels of acceptable LD between included SNPs,
min. number of individuals with data, and possible enrichment for "ancestry-informativeness" using min. # individuals 
and allele frequency differencess in allopatric maize and mex reference panels

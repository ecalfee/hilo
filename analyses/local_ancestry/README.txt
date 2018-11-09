Local ancestry inference for each admixed individual is done using ancestry_hmm
1) Thin sites to moderately low LD (spacing>.001cM ~ 1kb) and high ancestry informativeness for local ancestry inference
1a) calcAlleleFreqPop.sh: get allele frequencies uisng ANGSD at all sites for reference maize and reference mexicana panel
1b) filterAlloMAFnInd.sh: uses R script filter_var_sites_min_allo.R to filter out sites not meeting minimum # allopatric individuals
1c) pruneFixedcM.sh: uses python script to create a new sites file with a set of SNPs in low LD 
with high diff. in allele freq between ref panels and meeting minimum # individuals with data in each ref. panel
(see estimate in plot_mapping_metrics.R) minInd 5 and 11 are similar cutoffs ~86% sites should pass for mex or maize
But looking at an actual file, minInd 4 for mexicana cuts out ~1/3 of the data, so I'll use that as a bare minimum and somewhat arbitrarily minInd 10 for maize

Now I prepare input files for ancestry_hmm
2) Get allele counts for allopatric maize based on called genotypes:
2a) getVCFforAllo.sh: calls genotypes and makes a VCF for allopatric maize. 
Filters out genotypes for ind's with too high or too low coverage. 
Too high is defined as > 63 reads, this cutoff is
estimated in plot_mapping_metrics.R using depth of coverage for regions as > mu+3sd for highest coverage maize
2b) alloVCF2AlleleCounts.sh: gets major/minor allele counts out from VCF

3) Get major & minor reads for each hilo individual, including admixed and ref mexicana individuals
filters out data for individuals at sites where they have excess coverage
which is estimated as > 3sd from mean for the highest coverage hilo individual - over 8 reads  
countReadsACGT.sh: counts raw ACGT bases
countReadsMajorMinor.sh: associates ACGT with major/minor allele from sites file

4) Get allele counts for allopatric mexicana based on 1 sampled read per individual per site passing filtering:
4a) makeAlloCountsAncHMMInput.sh

5) Make an ancestry_hmm input file for each population
makeAncHMMInput.sh

6) Run ancestry_hmm with or without bootstrapping for time of admixture estimate
runAncHMM_noboot.sh
runAncHMM.sh

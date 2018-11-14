PCA and NGSadmix (= global ancestry) are used to visualize population structure
based on individual genotype likelihoods at a set of unlinked SNPs
1) Thin SNPs to low LD (<.01cM ~ 10kb) -- figure out how to do this by chromosome
pruneFixedcM.sh
variant sites results are stored in /thinnedPCA/
2) Calculate genotype likelihood at thinned SNPs using ANGSD
calcGL4Pruned.sh
3) Concatenate GL across all thinned SNPs & all chromosomes for input to NGSadmix
4) Run NGSadmix for global ancestry using K=2,3,4 genetic clusters
runNGSadmix.sh
5) Create PCA using PCAngsd and GL from all thinned SNPs
runPCAngsd.sh
6) Visualize results in R
7) Summarize population mean global ancestry for use in local ancestry calling,
omitting individuals with less than 0.05x mean coverage

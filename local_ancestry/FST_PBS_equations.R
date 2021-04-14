# Equations for calculating population differentiation:
# FST (2 populations) and Population Branch Statistic (PBS - 3 populations)

# Hudson 1992 pairwise Fst estimator from Bhatia 2013:
# (I implement equation 10 from the methods)
hudson_Fst <- function(p1, p2, n_hap1, n_hap2, min_hap = 2){ # allele frequencies and number of haplotypes (alleles) sampled
  if (n_hap1 < min_hap | n_hap2 < min_hap){
    fst = NA
  } else{
    N = (p1 - p2)^2 - p1*(1 - p1)/(n_hap1 - 1) - p2*(1-p2)/(n_hap2 - 1)
    D = p1*(1 - p2) + (1 - p1)*p2
    fst = N/D
  }
  return(fst)
}

# Population Branch Statistic from Yi 2013 10.1126/science.1190371
# (I implemented eq from supplement pg 8) and 
# apply using Hudson's 1992 Fst estimator above)

# branch length from fst:
branch_length = function(fst){
  -log(1 - fst)
}

# this gives the branch length specific to population 1:
pop_branch_stat = function(p1, p2, p3, n_hap1, n_hap2, n_hap3){
  if (p2 == p3 & abs(p1 - p3) == 1){
    # normally fixed SNPs would return NAs because they're undefined by the PBS, but we want to keep these very ancestry informative SNPs, so this clause
    pbs = Inf # returns high pbs values for fixed SNPs where pop 1 is fixed for a different allele than pops 2 & 3
  }
  else {
    fst12 = hudson_Fst(p1, p2, n_hap1, n_hap2)
    fst13 = hudson_Fst(p1, p3, n_hap1, n_hap3)
    fst23 = hudson_Fst(p2, p3, n_hap2, n_hap3)
    pbs = (branch_length(fst12) + branch_length(fst13) - branch_length(fst23))/2
  }
  return(pbs)
}

# test:
#is.na(pop_branch_stat(p1 = NA, p2 = 0, p3 = 0, n_hap1 = 10, n_hap2 = 10, n_hap3 = 10))
#is.na(pop_branch_stat(p1 = NA, p2 = NA, p3 = NA, n_hap1 = 10, n_hap2 = 10, n_hap3 = 10))
#is.na(pop_branch_stat(p1 = 1, p2 = 0, p3 = NA, n_hap1 = 10, n_hap2 = 10, n_hap3 = 10))
#hudson_Fst(p1 = 0, p2 = 0, n_hap1 = 10, n_hap2 = 10)
#hudson_Fst(p1 = 0.1, p2 = 0.1, n_hap1 = 1000, n_hap2 = 1000)
#pop_branch_stat(p1 = .6, p2 = 0.1, p3 = 0.1, n_hap1 = 100, n_hap2 = 100, n_hap3 = 100)
#pop_branch_stat(p1 = .6, p2 = 0.1, p3 = 0.1, n_hap1 = 100, n_hap2 = 100, n_hap3 = 100)
#pop_branch_stat(p1 = 1, p2 = 0.0001, p3 = 0.0001, n_hap1 = 100, n_hap2 = 100, n_hap3 = 100)
#pop_branch_stat(p1 = 1, p2 = 0.000000001, p3 = 0.000000001, n_hap1 = 100, n_hap2 = 100, n_hap3 = 100)
# note: we want to keep fixed snp differences as high pbs sites, even though they're traditionally excluded from pbs (e.g. Yi 2013)
#pop_branch_stat(p1 = 1, p2 = 0, p3 = 0, n_hap1 = 100, n_hap2 = 100, n_hap3 = 100)
#pop_branch_stat(p1 = 0, p2 = 1, p3 = 1, n_hap1 = 100, n_hap2 = 100, n_hap3 = 100)
#is.na(pop_branch_stat(p1 = 1, p2 = 0, p3 = 1, n_hap1 = 100, n_hap2 = 100, n_hap3 = 100))
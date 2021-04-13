# Hudson 1992 pairwise Fst estimator from Bhatia 2013 
# (I implement equation 10 from the methods;
# this code assumes 2 haplotypes observed per (diploid) individual observed) 
hudson_Fst <- function(p1, p2, n_ind1, n_ind2, min_ind = 2){ # allele frequencies and number of alleles sampled
  if (n_ind1 < min_ind | n_ind2 < min_ind){
    fst = NA
  } else{
    N = (p1 - p2)^2 - p1*(1 - p1)/(n_ind1*2 - 1) - p2*(1-p2)/(n_ind2*2 - 1)
    D = p1*(1 - p2) + (1 - p1)*p2
    fst = N/D
  }
  return(fst)
}

# Population Branch Statistic from Yi 2013 10.1126/science.1190371
# (I implemented eq from supplement pg 8) and 
# will apply using Hudson's 1992 Fst estimator)
# this gives the branch length specific to population 1:
# branch length from fst
branch_length = function(fst){
  -log(1 - fst)
}
pop_branch_stat = function(p1, p2, p3, n_ind1, n_ind2, n_ind3){
  fst12 = hudson_Fst(p1, p2, n_ind1, n_ind2)
  fst13 = hudson_Fst(p1, p3, n_ind1, n_ind3)
  fst23 = hudson_Fst(p2, p3, n_ind2, n_ind3)
  pbs = (branch_length(fst12) + branch_length(fst13) - branch_length(fst23))/2
}

# note: I want to keep fixed snp differences as high pbs, even though they're traditionally excluded from pbs (e.g. Yi 2013)



# convert sites to tracts
sites <- read.table("../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM/chr1.var.sites",
                    stringsAsFactors = F, header = F)
start = floor(sites$V2 - diff(c(sites$V2[1]-1, sites$V2), lag = 1)/2)
end = floor(sites$V2 + diff(c(sites$V2, sites$V2[length(sites$V2)]), lag = 1)/2)

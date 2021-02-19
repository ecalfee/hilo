library(dplyr)
library(tidyr)
#setwd("~/Documents/gitErin/hilo")
# load variables from Snakefile
Ne = snakemake@params[["Ne"]]
# Ne = 10000
yesno = snakemake@params[["yesno"]]
# yesno = "yes"
prefix_all = snakemake@params[["prefix"]]
# prefix_all = "HILO_MAIZE55"
txt_out = snakemake@output[["txt"]] # txt output file
# txt_out = paste0("ZAnc/results/", prefix_all, "/Ne", Ne, "_", yesno, "Boot/flowering_time_genes_v4.plus20kb.overlap.summary_overlap_outliers.txt")

mex_maize = c("mexicana", "maize")
posneg = c("pos", "neg")
stats = c("meanAnc", "lmElev")
sigs = c("fdr05", "perc02", "p05")

# function to calculate p-values for overlap based on permutation test
calc_pval = function(zea, sign, stat, sig, prefix_all, Ne, yesno){
  results_file = paste0("ZAnc/results/", prefix_all, "/Ne", Ne, "_", yesno, 
                        "Boot/flowering_time_genes_v4.plus20kb.", zea, "_", sign,
                        "_", stat, "_outliers.", sig, ".counts")
  count = read.table(results_file, header = F, nrows = 1)[1,1]
  shuffle = read.table(results_file, header = F, skip = 1)$V1
  pval = sum(shuffle >= count)/length(shuffle) # how many random shuffles have >= overlapping genes?
  return(data.frame(zea = zea, sign = sign, stat = stat, sig = sig, gene_count = count, pval = pval,
                    stringsAsFactors = F))
}

all_combs = expand.grid(zea = mex_maize, sign = posneg, stat = stats, sig = sigs)

results = do.call(rbind,
                  apply(all_combs, 1, function(row) 
  calc_pval(zea = row[["zea"]], sign = row[["sign"]], stat = row[["stat"]], sig = row[["sig"]], 
                  prefix_all = prefix_all, Ne = Ne, yesno = yesno)))

write.table(results, txt_out, sep = "\t", col.names = T, row.names = F, quote = F)

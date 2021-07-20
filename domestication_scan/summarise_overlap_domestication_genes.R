library(dplyr)
library(tidyr)
#setwd("~/Documents/gitErin/hilo")
# load variables from Snakefile
Ne = snakemake@params[["Ne"]]
# Ne = 10000
yesno = snakemake@params[["yesno"]]
# yesno = "yes"
prefix = snakemake@params[["prefix"]]
# prefix = "HILO_MAIZE55_PARV50"
K = snakemake@params[["K"]]
# K = 3
txt_out = snakemake@output[["txt"]] # txt output file
# txt_out = paste0("domestication_scan/results/", prefix, "/K", K, "/Ne", Ne, "_", yesno, "Boot/domestication_genes_from_lit.plus20kb.overlap.summary_overlap_outliers.txt")

# function to calculate p-values for overlap based on permutation test
calc_pval = function(zea, sign, stat, sig, prefix, K, Ne, yesno){
  results_file = paste0("domestication_scan/results/", prefix, "/K", K, "/Ne", Ne, "_", yesno, 
                        "Boot/domestication_genes_from_lit.plus20kb.", zea, "_", sign,
                        "_", stat, "_outliers.", sig, ".counts")
  count = read.table(results_file, header = F, nrows = 1)[1,1]
  shuffle = read.table(results_file, header = F, skip = 1)$V1
  pval = sum(shuffle >= count)/length(shuffle) # how many random shuffles have >= overlapping genes?
  return(data.frame(zea = zea, sign = sign, stat = stat, sig = sig, gene_count = count, pval = pval,
                    stringsAsFactors = F))
}

all_combs = data.frame(zea = c("mexicana", "maize"), 
                       sign = c("neg", "pos"), 
                       stat = "maize_anc", 
                       sig = "perc05")

results = do.call(rbind,
                  apply(all_combs, 1, function(row) 
  calc_pval(zea = row[["zea"]], sign = row[["sign"]], stat = row[["stat"]], sig = row[["sig"]], 
                  prefix = prefix, K = K, Ne = Ne, yesno = yesno)))

write.table(results, txt_out, sep = "\t", col.names = T, row.names = F, quote = F)

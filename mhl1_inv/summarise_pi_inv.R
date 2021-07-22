#!/usr/bin/env Rscript
library(dplyr)
library(xtable)

# this script summarises pi within zea groups
# for the maize and mexicana allele clusters
# across the putative inversion

# load variables from Snakefile
colors_file = snakemake@params[["colors"]]
# colors_file = "colors.R"
PREFIX = snakemake@params[["PREFIX"]]
# PREFIX = "HILO_MAIZE55_PARV50"
K = snakemake@params[["K"]]
# K = 3
Ne = snakemake@params[["Ne"]]
# Ne = 10000
YESNO = snakemake@params[["YESNO"]]
# YESNO = "yes"
tex_out = snakemake@output[["tex"]]
# tex_out = paste0("../hilo_manuscript/tables/", PREFIX, "_K", K, "_Ne", Ne, "_", YESNO, "Boot_summary_pi_inv_mhl1.tex")
tbl_out = snakemake@output[["tbl"]]
# tbl_out = paste0("mhl1_inv/results/", PREFIX, "/K", K, "/Ne", Ne, "_", YESNO, "Boot/summary_pi_inv_mhl1.txt")

# load data
source(colors_file)
ancestries = c("maize", "maize", "mexicana", "mexicana")
zea = c("parv", "maize", "maize", "mexicana")
names = c("parviglumis", "maize", "maize", "mexicana")
groups = paste(names, "within", ancestries, "inversion cluster")

d = do.call(bind_rows, lapply(1:4, function(i)
  read.table(paste0("mhl1_inv/results/", PREFIX, "/K", K, "/Ne", Ne, "_", YESNO, "Boot/", 
                      ancestries[i], "_cluster/", zea[i], ".pi.inv.pestPG"), 
header = F, skip = 1, sep = "\t") %>%
  data.table::setnames(c("region", "Chr",	"WinCenter",	"tW",	"tP", "tF",	"tH",	"tL",	"Tajima",	"fuf",	"fud",	"fayh",	"zeng",	"nSites")) %>%
  dplyr::select("tW", "tP", "nSites") %>%
  dplyr::mutate(wattersons_theta = tW/nSites,
                pairwise_theta = tP/nSites,
                inv_cluster = ancestries[i],
                zea = zea[i],
                sample = groups[i]))) %>%
  dplyr::select(., sample, pairwise_theta, wattersons_theta) %>%
  arrange(-pairwise_theta)

write.table(d, file = tbl_out, col.names = T, row.names = F, quote = F, sep = "\t")

print(xtable(d,
             type = "latex",
             latex.environments = NULL,
             digits = 3),
      include.rownames = F,
      file = tex_out)


#!/usr/bin/env Rscript
# script that add parviglumis (PARV50) to metadata
library(dplyr)
library(tidyr)

# snakemake input files
ids_file = snakemake@input[["ids"]]
# ids_file = "samples/HILO_MAIZE55_PARV50_ids.list"
meta_in = snakemake@input[["meta"]]
# meta_in = "samples/HILO_MAIZE55_meta.RData"

# output file
txt_out = snakemake@output[["txt"]]
# txt_out = "samples/HILO_MAIZE55_PARV50_meta.txt"
rdata_out = snakemake@output[["rdata"]]
# rdata_out = "samples/HILO_MAIZE55_PARV50_meta.RData"


ids <- read.table(ids_file, header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("ID"))
load(meta_in) # loads 'meta' data frame
meta_hilo_maize55 <- meta
meta <- ids %>% # keep same order of as input ids
  left_join(., meta_hilo_maize55, by = "ID") %>%
  dplyr::mutate(is_parviglumis = substr(ID, 1, 4) == "PARV",
                popN = ifelse(is_parviglumis, 3000, popN),
                zea = ifelse(is_parviglumis, "parviglumis", zea), 
                GEOCTY = ifelse(is_parviglumis, "Mexico", GEOCTY),
                group = ifelse(is_parviglumis, "parviglumis", group),
                LOCALITY = ifelse(is_parviglumis, "Palmar Chico", LOCALITY),
                dataset = ifelse(is_parviglumis, "PARV50", dataset)) %>%
  dplyr::select(-is_parviglumis)
# table(ids$ID == meta$ID) # all TRUE because order is preserved

# write output
write.table(meta, txt_out, 
            col.names = T, row.names = F,
            sep = "\t", quote = F)
save(list = c("meta"), file = rdata_out)
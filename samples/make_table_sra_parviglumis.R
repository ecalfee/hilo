#!/usr/bin/env Rscript
# working directory is hilo/
# make table of population metadata for manuscript
library(dplyr)
library(tidyr)

# from RI lab - 50 sample IDs we included
parv50 <- read.table("data/parviglumis/SRA_parviglumis_50_IDs.txt",
                   stringsAsFactors = F, header = F)$V1

# from ncbi
sra <- read.table("data/parviglumis/SraRunTable.txt", sep = ",",
                            stringsAsFactors = F, header = T) %>%
  dplyr::select(Run, Isolate) %>%
  dplyr::filter(Isolate %in% parv50)

write.table(sra, file = "samples/parviglumis_50_SRA_IDs.csv",
            sep = ",", col.names = T, row.names = F, quote = F)
write.table(sra, file = "../hilo_manuscript/files/parviglumis_50_SRA_IDs.csv",
            sep = ",", col.names = T, row.names = F, quote = F)

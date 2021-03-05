#!/usr/bin/env Rscript
# working directory is hilo/
# make table of population metadata for manuscript
library(dplyr)
library(tidyr)

load("samples/HILO_MAIZE55_meta.RData")
pops <- meta %>%
  dplyr::filter(group != "allopatric_maize") %>%
  dplyr::select(zea, symp_allo, LOCALITY, GEOCTY, ELEVATION, LAT, LON, RI_ACCESSION) %>%
  dplyr::filter(!duplicated(.)) %>%
  dplyr::rename(subspecies = zea,
         location = LOCALITY,
         country = GEOCTY,
         elevation = ELEVATION,
         latitude = LAT,
         longitude = LON,
         accession = RI_ACCESSION) %>%
  arrange(symp_allo, subspecies, elevation)

write.table(pops, file = "samples/population_metadata.csv",
            sep = ",", col.names = T, row.names = F, quote = F)

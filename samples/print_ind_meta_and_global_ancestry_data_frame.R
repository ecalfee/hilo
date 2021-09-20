# script loads R data and creates a tidy csv file
# for others to identify population, global ancestry estimate (NGSAdmix),
# for each individual sample
# working directory is hilo/
load("global_ancestry/results/NGSAdmix/HILO_MAIZE55_PARV50/K3_alphas_by_ind.RData")
d <- d_admix2 %>%
  dplyr::select(., ID, popN, RI_ACCESSION,
                zea, 
                #symp_allo, group, 
                LOCALITY,
                #ELEVATION, LAT, LON, 
                reads_q30,
                est_coverage, 
                #dataset, 
                parviglumis,mexicana, maize) %>%
  dplyr::rename(parviglumis_ancestry = parviglumis,
                mexicana_ancestry = mexicana,
                maize_ancestry = maize) %>%
  dplyr::mutate(has_local_ancestry = (est_coverage >= 0.5 & !reference_sample),
                has_local_ancestry = ifelse(reference_sample, NA, has_local_ancestry))
write.table(d, "samples/ind_metadata_and_genomewide_ancestry.csv",
            col.names = T, row.names = F, sep = ",")

d_pop <- read.table("samples/population_metadata.csv",
                    header = T, sep = ",") %>%
  dplyr::mutate(sample_type = ifelse(symp_allo == "allopatric",
                              "reference",
                              "sympatric")) %>%
  dplyr::select(-symp_allo) %>%
  bind_rows(.,
            data.frame(subspecies = c("maize", "parviglumis"),
                       location = "Palmar Chico",
                       country = "Mexico",
                       elevation = c(983, 1008),
                       sample_type = "reference")) %>%
  dplyr::arrange(sample_type, elevation, subspecies)

write.table(d_pop, "samples/population_metadata_HILO_MAIZE55_PARV50.csv",
            col.names = T, row.names = F, sep = ",")

write.table(d_pop, "../hilo_manuscript/files/S1_Table.csv",
            col.names = T, row.names = F, sep = ",")

library(dplyr)

samples <- read.csv("seeds_to_plant_5.1.19_new_IDs.txt",
                      stringsAsFactors = F, header = T, sep = "\t") %>%
  mutate(family = substr(family, 1, 1)) # so E1 becomes just E
water_1 <- read.csv("Hilo_Erin_plate1_amnt_water.csv",
                    stringsAsFactors = F, header = T, sep = ",") %>%
  mutate(plate = "one")
water_a <- read.csv("Hilo_Erin_platealpha_amnt_water.csv",
                    stringsAsFactors = F, header = T, sep = ",") %>%
  mutate(plate = "alpha")
water <- bind_rows(water_1, water_a) %>%
  mutate(letter = substr(plate_well_code, 1, 1),
         number = as.integer(substr(plate_well_code, 2, 3))) %>%
  dplyr::select(plate, conc, letter, number)

lab_id <- read.csv("get_tissue_5.28.19.txt",
                   stringsAsFactors = F, header = T, sep = "\t") %>%
  mutate(family = substr(family, 1, 1)) %>%
  mutate(plate = ifelse(plate %in% c("X", "Y", "Z"), "one", plate)) 
lab_id_water <- dplyr::left_join(lab_id, water, by = c("plate", "letter", "number"))

combined <- full_join(samples, lab_id_water, 
                      by = c("RI_ACCESSION", "popN", "family")) %>%
  filter(is.na(Exclude.repeat)) %>%
  mutate(expected_coverage = 
           ifelse(!is.na(conc), 1.5, est_coverage))
# likely library fails:
combined %>%
  filter(!is.na(Library.fail) | conc <= 2) %>%
  #dplyr::select(Library.fail, conc) 
  group_by(popN, LOCALITY, zea, symp_allo) %>%
  summarise(tot_coverage = sum(expected_coverage),
            n_ind_include = sum(expected_coverage >= 0.5)) %>%
  View(.)
# likely library successes:
combined %>%
  filter(is.na(Library.fail) & expected_coverage >= 0.5) %>%
  filter(is.na(conc) | conc > 2) %>%
  group_by(popN, LOCALITY, zea, symp_allo) %>%
  summarise(tot_coverage = sum(expected_coverage),
            n_ind_include = sum(expected_coverage >= 0.5)) %>%
  View(.)

# expected totals after reseq:
combined %>%
  filter(is.na(Library.fail)) %>% # & expected_coverage >= 0.5) %>%
  filter(is.na(conc) | conc > 0.4) %>%
  filter(expected_coverage >= 0.5) %>% # note: meta data loaded from plot_NGSAdmix.R
  bind_rows(., mutate(meta[meta$ID %in% landraces$ID, ], expected_coverage = est_coverage)) %>%
  filter(LOCALITY != "Puerta Encantada") %>%
  #group_by(LOCALITY, zea, symp_allo) %>%
  group_by(zea, symp_allo) %>%
  summarise(tot_coverage = sum(expected_coverage),
            n_ind_over0.5x = sum(expected_coverage >= 0.5),
            n_ind_over1x = sum(expected_coverage >= 1),
            mean_coverage = mean(expected_coverage),
            max_coverage = max(expected_coverage)) %>%
  arrange(symp_allo, zea) %>%
  View(.)


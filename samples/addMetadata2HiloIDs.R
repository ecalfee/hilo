# made new script that links metadata to hilo ID using updated ID-population link from Anne 8.8.18
library(dplyr)
library(tidyr)
library(rgbif) # for elevation data
source("../colors.R")

# ID and population labels:
pools <- do.call(rbind,
                 lapply(1:8, function(p) read.csv(paste0("../samples/pool", p, ".csv"),
                                                  stringsAsFactors = F)[,1:7] %>%
                          mutate(., lane = p) %>%
                          filter(!(Library.name=="")))) %>%
  dplyr::rename(ID = Library.name) %>%
  separate(data = ., col = sample_name, sep = "_", c("popN", "family"), extra = "merge") %>%
  mutate(popN = as.integer(popN)) %>%
  dplyr::select(c("ID", "Adapter", "lane", "plate_number", "popN", "family"))
pools$plate_number[pools$plate_number=="plate 3+4"] <- "plate3+4" # plate label inconsistency
# fix 61-66 label switch:
old61 <- pools$popN[pools$ID=="HILO61"]
old66 <- pools$popN[pools$ID=="HILO66"]
pools$popN[pools$ID=="HILO61"] <- old66
pools$popN[pools$ID=="HILO66"] <- old61

meta <- pools %>%
  mutate(., n = as.integer(substr(ID, 5, 100))) %>%
  dplyr::select(., c("ID", "n", "Adapter", "lane", "plate_number", "popN", "family")) %>%
  mutate(., popN = as.numeric(popN)) %>%
  mutate(., zea = ifelse(popN >= 100, "maize", "mexicana")) %>%
  mutate(., symp_allo = ifelse(popN %in% c(20, 22, 33), "allopatric", "sympatric")) %>%
  dplyr::arrange(ID)

# add sequencing round
seq_round <- read.table("../data/HILO_raw_reads/merged_all.list",
                        header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("ID", "seq_lane", "library")) 
seq_round2 <- seq_round %>%
  bind_cols(., filter(pools, !(ID == "HILO80" & lane == 2))) %>%
  dplyr::rename(p7 = Adapter,
                fam_ind = family) %>%
  dplyr::mutate(instrument_model = "Illumina HiSeq 4000", # add sequencing machine
                sample_name = ID,
                incl_excl = ifelse((ID == 80 | plate_number == "plate2"), 
                                   "excl_label_error", NA), # mark excluded samples
                fam_ind = stringr::str_trim(fam_ind, side = "both"), # trim white space
                family = stringr::str_sub(fam_ind, 1, 1),
                pop_fam = paste(popN, family, sep = "_"),# I arbitrarily assign x and xi to these individuals -- tissue was not saved
                ind = ifelse(duplicated(pop_fam) & !(duplicated(ID)), "x", "xi") # will be NA if no number
                ) # now add RIMME/RIMMA.
# basic checks:
seq_round2 %>%
  group_by(pop_fam) %>%
  summarise(ids = length(unique(ID)),
            inds = length(unique(ind))) %>%
  filter(ids != 1 | inds != 1) # check -- good, some families are repeated but have multiple individual ids
seq_round2 %>%
  group_by(ID) %>%
  summarise(inds = length(unique(paste(pop_fam, ind)))) %>%
  filter(inds != 1) # check -- good, each ID has only 1 unique pop_fam_ind identifier
table(seq_round$ID == seq_round$ID1) # good
unique(paste0(seq_round2$seq_lane, "_", seq_round2$lane)) # good lanes 1-5 were March2018, 7-8 were Jan2019

# add in latest novaseq sequencing data
april2020 <- read.table("../data/HILO_raw_reads/April2020_all.list", stringsAsFactor = F,
                        header = F) %>%
  data.table::setnames(c("ID", "seq_lane", "library")) %>%
  left_join(.,
                     read.table("../samples/april2020_meta_from_Taylor.txt", 
                                header = T, stringsAsFactors = F),
                     by = "ID") %>%
  mutate(rimme = stringr::str_sub(pop_family, 1, 9),
         family = stringr::str_sub(pop_family, -1, -1),
         popN = as.integer(stringr::str_sub(pop_family, -4, -2)),
         lane = ifelse(seq_lane == "April2020_1", 9, 10),
         plate = ifelse(seq_lane == "April2020_1", "regrow_1", "regrow_alpha"),
         instrument_model = "Illumina HiSeq 6000")
backup <- read.table("get_tissue_5.28.19.txt", header = T, stringsAsFactors = F,
                     sep = "\t") %>%
  dplyr::rename(tissue_collection_plate = plate) %>%
  mutate(well = paste0(letter, number)) %>% # just check that it matches. then make 1-> i and 2->ii
  mutate(plate = ifelse(tissue_collection_plate == "alpha", 
                        "regrow_alpha", "regrow_1"))
novaseq <- backup %>%
  dplyr::mutate(family = stringr::str_sub(family, 1, 1)) %>%
  dplyr::select(RI_ACCESSION, plate, well, popN, family, tissue_collection_plate,
                ind2seq, regrow_batch) %>%
  full_join(., april2020, by = c("well", "plate", "RI_ACCESSION"="rimme", "popN", "family")) %>%
  dplyr::filter(!is.na(ID)) %>%
  mutate(ind = ifelse(ind2seq == 1, "i", "ii"))

# all plates
all_libraries <- seq_round2 %>%
  dplyr::rename(plate = plate_number) %>%
  mutate(RI_ACCESSION = paste0(ifelse(popN < 100, "RIMME", "RIMMA"),
                               stringr::str_pad(popN, width = 4, "left", pad = "0"))) %>%
  bind_rows(., novaseq) %>%
  mutate(platform = "ILLUMINA",
         library_ID = paste0(ID, "_L", stringr::str_pad(lane, width = 3, "left", pad = "0"))) %>%
  dplyr::select(ID, RI_ACCESSION, popN, family, ind, incl_excl, library_ID, library, lane, seq_lane, plate, p7, well, regrow_batch, platform, instrument_model)

# basic checks:
all_libraries %>%
  mutate(pop_fam = paste(popN, family, sep = "_")) %>%
  group_by(pop_fam) %>%
  summarise(ids = length(unique(ID)),
            inds = length(unique(ind))) %>%
  filter(ids != inds)# check -- good, some families are repeated but have multiple individual ids


# add in germplasm data from JRI for all projects ("riplasm").

# write out full file
# get list of all the file names from farm
# 
# check how I did file merging previously from multiple bams -- 
# do I need to do this or should I just ditch the very low coverage one? (or both)

# RIMMA is hilo maize and RIMME is hilo mex. This just gives me metadata for each population..(e.g. lat/long)
riplasm <- read.csv("../data/riplasm/riplasm.csv", stringsAsFactors = F, header = T) %>%
  mutate(., prefix = substr(RI_ACCESSION, 1, 5)) %>%
  mutate(., ind = as.integer(substr(RI_ACCESSION, 6, 100))) %>%
  dplyr::filter(., prefix %in% c("RIMMA", "RIMME"))

# sticking it all together
hilo <- cbind(meta, do.call(rbind,
                             apply(meta, 1, function(row) riplasm[riplasm$RI_ACCESSION ==  paste0(
                               ifelse(row["zea"] == "mexicana", "RIMME", "RIMMA"),
                               sprintf("%04d", as.integer(row["popN"]))),
                               c("RI_ACCESSION", "GEOCTY", "LOCALITY")]))) %>%
  dplyr::arrange(., n)

#unique(hilo$LOCALITY) - some inconsistency in naming leads to false repeats
hilo$LOCALITY[hilo$LOCALITY == "Puerta Encantada (Amatlan)"] <- "Puerta Encantada"
hilo$LOCALITY[hilo$LOCALITY == "Jicaltepec (San Pablo Autopan)"] <- "Jicaltepec"
hilo$LOCALITY[hilo$LOCALITY == "Amatlan, Morelos"] <- "Amatlan"
hilo$LOCALITY[hilo$LOCALITY == "Puruandiro(Victor)"] <- "Puruandiro"

# write file with all hilo individuals
dplyr::select(hilo, c("ID", "n", "popN", "family", "zea", "symp_allo", "RI_ACCESSION", "GEOCTY", "LOCALITY")) %>%
  unique(.) %>%
  write.table(., "hilo_meta.txt", row.names = F, quote = F, col.names = T, sep = "\t")

dplyr::select(hilo, c("ID", "Adapter", "lane", "plate_number")) %>%
                write.table(., "hilo_plates.txt", row.names = F, quote = F, col.names = T, sep = "\t")

# write ID file for samples included in 'pass2'
# pass2 excludes PCR plate #2 for label switch/permutation and HILO80 based on lab note of mixup or contamination of that well
hilo$pass2 <- !(hilo$ID == "HILO80" | hilo$plate_number == "plate2")
# write pass2 file
write.table(unique(hilo$ID[hilo$pass2 == T]), "pass2_IDs.list", row.names = F, col.names = F, quote = F)

# make population files for pass2 included samples
pass2.ind <- hilo %>% # only include individuals once, even if sequenced twice
  filter(., pass2 == T) %>%
  dplyr::select(., c("ID", "popN", "family", "zea", "symp_allo", "RI_ACCESSION", "GEOCTY", "LOCALITY")) %>%
  unique(.)
sum(duplicated(pass2.ind$ID)) # good no repeat IDs
# make population files - sympatric maize and mexicana
for (z in c("maize", "mexicana")){
  pass2.ind %>%
    filter(., zea == z) %>%
    filter(., symp_allo == "sympatric") %>%
    dplyr::select(., "ID") %>%
    write.table(.,
                paste0("pass2_pops/symp.", z, "_IDs.list"),
                row.names = F, col.names = F, quote = F)

}
# allopatric mexicana
pass2.ind %>%
  filter(., zea == "mexicana") %>%
  filter(., symp_allo == "allopatric") %>%
  dplyr::select(., "ID") %>%
  write.table(.,
              paste0("pass2_pops/allo.mexicana_IDs.list"),
              row.names = F, col.names = F, quote = F)

# all other populations
for (N in unique(pass2.ind$popN)){
  pass2.ind %>%
    filter(., popN == N) %>%
    dplyr::select(., "ID") %>%
    write.table(.,
                paste0("pass2_pops/pop", N, "_IDs.list"),
                row.names = F, col.names = F, quote = F)

}


# make metadata file for 4 lowland maize populations from mexico populations:
maize4low <- data.frame(ID = read.table("MAIZE4LOW_IDs.list", header = F, stringsAsFactors = F)$V1,
           popN = c(0,0,0,1,1,1,1,2,2,2,2,3,3,3,3),
           GEOCTY = "Mexico",
           zea = "maize",
           symp_allo = "allopatric") %>%
  dplyr::mutate(LOCALITY = paste("mex4low", popN, sep = "_")) %>%
  dplyr::mutate(n = substr(ID, 10, 12))
# write file
write.table(maize4low, "MAIZE4LOW_meta.txt", row.names = F, quote = F, col.names = T, sep = "\t")

# make metadata file for allopatric maize landraces:
landraces <- data.frame(ID = read.table("alloMAIZE_IDs.list", header = F, stringsAsFactors = F)$V1,
                        zea = "maize",
                        symp_allo = "allopatric",
                        n = 1:15) %>%
  left_join(., dplyr::select(riplasm, c(GEOCTY, RI_ACCESSION)) %>% mutate(ID = RI_ACCESSION), by = "ID") %>%
  mutate(., popN = as.integer(substr(ID, 6, 9)))
# filling in missing/corrections for the landrace IDs
landraces$GEOCTY[landraces$ID == "RIMMA0703"] <- "Mexico" # clear from listed locality, Mexico Federal District
landraces$GEOCTY[landraces$ID == "RIMMA0625"] <- riplasm[riplasm$RI_ACCESSION=="RIMMA0438", "GEOCTY"] # label switch ID-ed by Li
landraces$LOCALITY <- ifelse(landraces$GEOCTY %in% c("Mexico", "Guatemala"), "MexLow",
                             ifelse(landraces$GEOCTY %in% c("Peru", "Ecuador"), "Andes",
                                    ifelse(landraces$GEOCTY=="Colombia", "SALow", NA)))

# write file
write.table(landraces, "alloMAIZE_meta.txt", row.names = F, quote = F, col.names = T, sep = "\t")


# get location information and elevation
# add elevation data for populations/locations also used in Hufford 2013
elevHufford2013Supp = data.frame(LOCALITY = c("Jicaltepec", "Puruandiro", "Nabogame",
                                              "Penjamillo", "Tlapala", "Cocotilan", "San Pedro",
                                              "Opopeo", "Amatlan", "Xochimilco", "Tenango del Aire",
                                              "Puerta Encantada", "Santa Clara", "Amecameca",
                                              "Malinalco", "Ixtlan", "El Porvenir"),
                                 ELEVATION = c(NA, 1915, 2020,
                                               NA, NA, NA, 2459,
                                               2213, 1658, 2237, 2609,
                                               1658, 2173, NA,
                                               NA, 1547, 2094), stringsAsFactors = F) %>%
  filter(., !is.na(ELEVATION))
# accessed google api through web on Aug 29 2018 to retrieve missing elevation
elevGoogleAPI2018 = data.frame(LOCALITY = c("Amecameca", "Jicaltepec",
                                            "Penjamillo", "Cocotilan",
                                            "Malinalco", "Tlapala"),
                               ELEVATION = c(2467, 2587,
                                             1705, 2269,
                                             1887, 2272), stringsAsFactors = F)
gps <- unique(hilo[ , c("popN", "zea", "symp_allo", "RI_ACCESSION", "GEOCTY", "LOCALITY")]) %>%
  left_join(., riplasm[ , c("RI_ACCESSION", "LAT", "LON")], by = "RI_ACCESSION") %>%
  left_join(., rbind(elevHufford2013Supp, elevGoogleAPI2018),
            by = "LOCALITY") %>%
  mutate(., LON = -1*LON) # original data has degree decimal's west as positive not negative numbers

apikey <- "AIzaSyBLqWTiXQEG4m9Vu1oVmdvHke5sMw9gp6Y"
# very limited number of API request per day -- first look for missing values
elev1 = unique(gps[is.na(gps$ELEVATION), c("LOCALITY", "LAT", "LON")]) %>%
  mutate(., decimalLatitude = LAT) %>%
  mutate(., decimalLongitude = LON) %>%
  elevation(input = .,
            key = apiKey)
# then confirm existing values
elev2 = unique(gps[!is.na(gps$ELEVATION), c("LOCALITY", "LAT", "LON")]) %>%
  mutate(., decimalLatitude = LAT) %>%
  mutate(., decimalLongitude = LON) %>%
  elevation(input = .,
                 key = apiKey)

write.table(gps, "gps_and_elevation_for_sample_sites.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")

# add metrics data for pass2 coverage:
metrics <- read.table("../filtered_bams/metrics/hilo_alloMAIZE_MAIZE4LOW.flagstat.total", stringsAsFactors = F)
colnames(metrics) <- c("ID", "total_reads_pass")
# size of reference genome reads are mapped to
ref_genome_size <- sum(read.table("../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.fai",
                                    stringsAsFactors = F)$V2) # V2 is the size of each chromosome mapped to,
# including parts I won't analyze on the Pt and Mt and scaffolds not assigned to chromosomes (but reads would pass Q filters there too)
metrics$est_coverage = round(metrics$total_reads_pass*150/ref_genome_size, 4)

length(unique(hilo$ID))
length(unique(hilo$Adapter))
dim(unique(hilo[, c("ID", "Adapter")]))
hilo[!duplicated(hilo[, c("ID", "Adapter")]) & duplicated(hilo$ID), c("ID", "Adapter")]
# the only sample with more than one adapter is HILO80; which is excluded

dim(unique(hilo[, c("ID", "Adapter")]))

target = 8
status <- pass2.ind %>%
  left_join(., metrics, by = "ID") %>%
  group_by(popN, LOCALITY, RI_ACCESSION, symp_allo, zea) %>%
  summarise(over0.5x = sum(est_coverage >= 0.5),
            reseq0.25to0.5x = sum(est_coverage < 0.5 & est_coverage >= 0.25),
            under0.25x = sum(est_coverage < 0.25)) %>%
  mutate(need_to_plant = target - over0.5x - reseq0.25to0.5x) %>%
  mutate(need_to_plant = ifelse(need_to_plant < 0, 0, need_to_plant))

#cov_below_0.25x <- pass2.ind %>%
pass2.ind %>%
  left_join(., metrics, by = "ID") %>%
  #filter(est_coverage < 0.25) %>%
  ggplot(., aes(x = est_coverage, fill = paste(zea, symp_allo, sep = "_"))) +
  geom_histogram()
hist(cov_below_0.25x$est_coverage, breaks = 20)
pass2.ind %>%
  left_join(., metrics, by = "ID") %>%
  filter(est_coverage < 0.25) %>%
  ggplot(., aes(x = est_coverage, fill = paste(zea, symp_allo, sep = "_"))) +
  geom_histogram()
pass2.ind %>%
  left_join(., metrics, by = "ID") %>%
  filter(., est_coverage >= 0.25 & est_coverage < .5) %>%
  ggplot(., aes(x = est_coverage*2, fill = paste(zea, symp_allo, sep = "_"))) +
  geom_histogram()
pass2.ind %>%
  left_join(., metrics, by = "ID") %>%
  #filter(., est_coverage < .5) %>%
  filter(., symp_allo == "allopatric") %>%
  ggplot(., aes(x = est_coverage, fill = LOCALITY)) +
  geom_histogram()

pass2.ind %>%
  left_join(., unique(hilo[ , c("ID", "Adapter")]), by = "ID") %>%
  #filter(symp_allo != "sympatric") %>%
  #filter(popN != 20) %>%
  left_join(., metrics, by = "ID") %>%
  filter(., est_coverage >= 0.2 & est_coverage < .5) %>%
  #group_by(lane) %>%
  group_by(Adapter) %>%
  summarise(n()) %>%
  View(.) # great, no more than 2 uses per adapter so I can repool and split across 2 lanes

sum(status$need_to_plant) + 5
sum(status$reseq0.25to0.5x)
(105+68)/4
(78+68)/4
View(status)
4*1300+500

status %>%
  filter(symp_allo == "allopatric") %>%
  View(.) # why are there > 10 for one pop?


# Coverage estimate to run 96 samples on a lane of NovaSeq6000:
n_pass2_sequencing <- nrow(pass2.ind) + sum(duplicated(hilo %>% # only include individuals once, even if sequenced twice
                                   filter(., pass2 == T) %>%
                                 dplyr::select(., ID)))
pass2.ind %>% left_join(., metrics, by = "ID") %>%
  filter(., !(zea == "maize" & symp_allo == "allopatric")) %>%
  dplyr::select(., est_coverage) %>%
  sum()/n_pass2_sequencing # .47-.55 mean coverage for 254 individuals; 40 per lane

sapply(c(.47, .55), function(x) c(2.5, 3)*10^9/(340*10^6)*x*40/96)
# so I can expect between 1.4 and 2x mean coverage per re-sequenced individual

# make a re-sequencing priority list:
filter(pass2.ind, symp_allo == "sympatric") %>%
  left_join(., metrics, by = "ID") %>%
  filter(est_coverage >= .5) %>%
  summarise(n())
length(unique(pass2.ind[pass2.ind$symp_allo == "sympatric", "RI_ACCESSION"]))*8
filter(pass2.ind) %>%
  left_join(., metrics, by = "ID") %>%
  filter(est_coverage >= .5) %>%
  ggplot(aes(x = LOCALITY, fill = paste(zea, symp_allo, sep = "_"))) +
  geom_bar() +
  facet_wrap(~zea)
hilo %>%
  dplyr::select(., c("ID", "zea", "symp_allo", "family", "popN", "LOCALITY", "pass2")) %>%
  unique(.) %>%
  ggplot(aes(x = LOCALITY, fill = paste(zea, symp_allo, sep = "_"))) +
  geom_bar() +
  facet_wrap(~zea)

reseq <- hilo %>%
  left_join(., metrics, by = "ID") %>%
  mutate(est_coverage = ifelse(ID == "HILO80" | plate_number == "plate2", 0, est_coverage)) %>%
  dplyr::select("RI_ACCESSION", "family", "zea", "symp_allo", "ID", "LOCALITY", "plate_number", "est_coverage") %>%
  unique()

flagstat <- read.table("../filtered_bams/metrics/flagstat/multiqc/multiqc_data/multiqc_samtools_flagstat.txt",
           header = T, stringsAsFactors = F, sep = "\t") %>%
  tidyr::separate(Sample, into = c("metrics", "flagstat", "seq_lane", "ID"),
                  sep = " \\| ") %>%
  dplyr::select(c("ID", "seq_lane", "total_passed")) %>%
  dplyr::rename(reads_q30 = total_passed)
  # divide sample into multiple columns, add population info, sort etc., ID exclusions
hilo2 <- left_join(all_libraries, flagstat, by = c("ID", "seq_lane")) %>%
  mutate(est_coverage = reads_q30*150/(ref_genome_size)) %>% # ~2.1*10^9 is the genome seq for the roughly 2.4Gb genome. this will be a slight underestimate since not all of the genome is mappable
  dplyr::filter(is.na(incl_excl)) %>%
  left_join(., dplyr::group_by(., ID) %>%
              dplyr::summarise(n_seq = n(), 
                               est_coverage_combined = sum(est_coverage)), # if I merge reads
            by = "ID") %>%
  left_join(., dplyr::group_by(., RI_ACCESSION, popN, family) %>%
              dplyr::summarise(best_ind = ind[which.max(est_coverage_combined)]),
            by = c("RI_ACCESSION", "popN", "family")) %>%
  arrange(desc(est_coverage_combined)) %>%
  mutate(best = (ind == best_ind),
         too_low = est_coverage_combined < .05,
         keep = best & !too_low)
hilo2 %>%
  filter(!duplicated(ID)) %>% # don't count samples twice
  dplyr::select(best, too_low) %>%
  table(.)

hilo2 %>%
  left_join(., gps) %>%
  filter(!duplicated(ID)) %>%
  filter(keep) %>%
  dplyr::group_by(popN, zea, symp_allo, LOCALITY, ELEVATION, RI_ACCESSION) %>%
  dplyr::summarise(pop_n = n(),
                   total_pop_depth = sum(est_coverage_combined)) %>%
  View(.)
hilo2 %>%
  left_join(., gps) %>%
  filter(!duplicated(ID)) %>%
  filter(keep) %>%
  filter(est_coverage_combined >=0.5) %>%
  filter(symp_allo == "sympatric") %>%
  dplyr::group_by(popN, zea, symp_allo, LOCALITY, ELEVATION, RI_ACCESSION) %>%
  dplyr::summarise(pop_n = n(),
                   total_pop_depth = sum(est_coverage_combined)) %>%
  dplyr::ungroup() %>%
  dplyr::summarise(min_n = min(pop_n),
                   max_n = max(pop_n),
                   tot_n = sum(pop_n))
hilo2 %>%
  left_join(., gps) %>%
  filter(!duplicated(ID)) %>%
  filter(keep) %>%
  filter(symp_allo == "allopatric") %>%
  summarise(n = n(),
            total_depth = sum(est_coverage_combined))


p_seq_counts <- hilo2 %>%
  left_join(., gps, by = c("RI_ACCESSION", "popN")) %>%
  filter(!duplicated(ID)) %>%
  mutate(coverage = cut(est_coverage_combined, 
                        breaks = c(0, 0.05, 0.5, 10)
                        )) %>%
  filter(keep) %>%
  ggplot(., aes(alpha = coverage,
                fill = zea,
                x = LOCALITY)) +
  geom_bar() +
  labs(fill = "Zea", x = "Population", y = "# Individuals") +
  scale_alpha_discrete(range = c(.6, 1),
                       name = "Coverage",
                       labels = c(expression("x "<" 0.5"), expression("x ">=" 0.5"))) +
  scale_fill_manual(values = col_maize_mex_parv) +
  facet_grid(zea~.) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))#
        #axis.title.x = "Population")
ggsave("../../hilo_manuscript/figures/p_seq_counts.png",
       plot = p_seq_counts,
       height = 4, width = 5.4, units = "in",
       device = "png")
hilo2 %>%
  filter(keep) %>%
  dplyr::summarise(min = min(est_coverage_combined),
    max = max(est_coverage_combined))

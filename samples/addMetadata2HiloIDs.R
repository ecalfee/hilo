# made new script that links metadata to hilo ID using updated ID-population link from Anne 8.8.18
library(dplyr)
library(tidyr)
#library(rgbif) # for elevation data
source("../colors.R")

# ID and population labels:
pools <- do.call(rbind,
                 lapply(1:8, function(p) read.csv(paste0("../samples/pools/pool", p, ".csv"),
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
# seq_round2 %>%
#   group_by(pop_fam) %>%
#   summarise(ids = length(unique(ID)),
#             inds = length(unique(ind))) %>%
#   filter(ids != 1 | inds != 1) # check -- good, some families are repeated but have multiple individual ids
# seq_round2 %>%
#   group_by(ID) %>%
#   summarise(inds = length(unique(paste(pop_fam, ind)))) %>%
#   filter(inds != 1) %>%
#   nrow(.) == 0 # check -- good, each ID has only 1 unique pop_fam_ind identifier
# unique(paste0(seq_round2$seq_lane, "_", seq_round2$lane)) # good lanes 1-5 were March2018, 7-8 were Jan2019

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
backup <- read.table("../samples/get_tissue_5.28.19.txt", header = T, stringsAsFactors = F,
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
  filter(ids != inds) %>%
  nrow(.) == 0 # check -- good, some families are repeated but have multiple individual ids

# add in germplasm data from JRI for all projects ("riplasm").
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

write.table(gps, "../samples/gps_and_elevation_for_sample_sites.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")

# get estimated coverage across maize genome
# size of reference genome reads are mapped to
ref_genome_size <- sum(read.table("../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.fai",
                                    stringsAsFactors = F)$V2) # V2 is the size of each chromosome mapped to,

# samtools flagstat read counts mapping quality > 30
flagstat <- read.table("../filtered_bams/metrics/flagstat/multiqc/multiqc_data/multiqc_samtools_flagstat.txt",
           header = T, stringsAsFactors = F, sep = "\t") %>%
  tidyr::separate(Sample, into = c("metrics", "flagstat", "seq_lane", "ID"),
                  sep = " \\| ") %>% # divide sample into multiple columns
  dplyr::select(c("ID", "seq_lane", "total_passed")) %>%
  dplyr::rename(reads_q30 = total_passed)
# add population info, sort etc., ID exclusions
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
         too_low = est_coverage_combined < .05, # exclude if < 0.05x coverage
         keep = best & !too_low) %>%
  left_join(., gps, by = c("popN", "RI_ACCESSION")) # add pop information
# how many samples are excluded vs. successfully resequenced with a replacement?
# hilo2 %>%
#   filter(!duplicated(ID)) %>% # don't count samples twice
#   dplyr::select(best, too_low) %>%
#   table(.)
# # total depth for allopatric maize
# hilo2 %>%
#   left_join(., gps) %>%
#   filter(!duplicated(ID)) %>%
#   filter(keep) %>%
#   filter(symp_allo == "allopatric") %>%
#   summarise(n = n(),
#             total_depth = sum(est_coverage_combined))

# included samples:
included <- hilo2 %>%
  filter(keep) %>%
  group_by(ID, popN, family, ind, RI_ACCESSION, zea, symp_allo, GEOCTY, LOCALITY, ELEVATION, LAT, LON) %>%
  summarise(reads_q30 = sum(reads_q30),
            est_coverage = sum(est_coverage)) %>%
  arrange(as.numeric(substr(ID, 5, 8))) # arrange by hilo ID
write.table(included, "../samples/HILO_meta.txt", col.names = T, row.names = F,
            sep = "\t", quote = F)

# write out included IDs and bams
hilo_ids = included$ID
write.table(hilo_ids, "../samples/HILO_ids.list",
            col.names = F, row.names = F, quote = F, sep= "\t")
hilo_bams = paste0("filtered_bams/merged_bams/", hilo_ids, ".sort.dedup.bam")
write.table(hilo_bams, 
            "../samples/HILO_bams.list",
            col.names = F, row.names = F, quote = F, sep= "\t")

# write out maize55 allopatric maize ids & bams (from Palmar Chico):
maize55_ids = paste0("MAIZE", 1:55)
write.table(maize55_ids, "../samples/MAIZE55_ids.list",
            col.names = F, row.names = F, quote = F, sep= "\t")
maize55_bams = paste0("filtered_bams/results/Maize55/", maize55_ids, ".sort.dedup.bam")
write.table(maize55_bams, "../samples/MAIZE55_bams.list",
            col.names = F, row.names = F, quote = F, sep= "\t")

# together
write.table(c(hilo_ids, maize55_ids), "../samples/HILO_MAIZE55_ids.list",
            col.names = F, row.names = F, quote = F, sep= "\t")
write.table(c(hilo_bams, maize55_bams), "../samples/HILO_MAIZE55_bams.list",
            col.names = F, row.names = F, quote = F, sep= "\t")


# write out fastq files & sequencing meta data for included samples:
fastq_info2 <- hilo2 %>%
  filter(keep) %>%
  dplyr::mutate(fastq_prefix = paste0("data/HILO_raw_reads/", seq_lane, "/", ID),
                fastq1 = paste0(fastq_prefix, "_1.fq.gz"),
                fastq2 = paste0(fastq_prefix, "_2.fq.gz"),
                pop_fam_ind = paste(RI_ACCESSION, family, ind, sep = "_")) %>%
  dplyr::select(ID, RI_ACCESSION, popN, family, ind, lane, 
                pop_fam_ind, library_ID, 
                platform, instrument_model, fastq1, fastq2) %>%
  arrange(as.numeric(substr(ID, 5, 8)), lane) # arrange by hilo ID
write.table(fastq_info2, "../samples/HILO_fastq_metadata.txt",
            col.names = T, row.names = F, sep = "\t")

p_seq_counts <- included %>%
  mutate(coverage = cut(est_coverage, 
                        breaks = c(0, 0.05, 0.5, 10)
                        )) %>%
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
#p_seq_counts
ggsave("../../hilo_manuscript/figures/p_seq_counts.png",
       plot = p_seq_counts,
       height = 4, width = 5.4, units = "in",
       device = "png")


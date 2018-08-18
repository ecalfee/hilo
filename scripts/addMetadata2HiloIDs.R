# made new script that links metadata to hilo ID using updated ID-population link from Anne 8.8.18
library(dplyr)
library(tidyr)
# I STILL NEED TO SWAP HILO66 AND HILO61 DUE TO AN ERROR THAT'S DOCUMENTED IN FIELD NOTES

# ID and population labels:

# meta 1 gives accurate metadata for the first 142 samples and came from Anne 
# over slack with the hiloID-to-popN link
meta1 <- read.csv("../data/pre_label_fix/HILO_samples.csv", stringsAsFactors = F, sep = ",") %>%
  select(., c("Library.name", "sample_name")) %>%
  separate(data = ., col = sample_name, sep = "_", c("popN", "family"), extra = "merge") %>%
  rename(., ID = Library.name) %>%
  mutate(., n = as.integer(substr(ID, 5, 100))) %>%
  filter(., n <= 142) # accurate pop data for ID's above 142 come from a second and third file
# fix population-level swap documented by Dan between HILO66 and HILO61
meta1_fixed = meta1
# separate population identifiers from incorrect ID
meta1_fixed[meta1_fixed$ID == "HILO66", c("popN", "family")] <- meta1[meta1$ID == "HILO61", c("popN", "family")]
meta1_fixed[meta1_fixed$ID == "HILO61", c("popN", "family")] <- meta1[meta1$ID == "HILO66", c("popN", "family")]

meta2 <- read.csv("../data/pre_label_fix/new_hiloID_link_from_Anne_8.8.18.txt", stringsAsFactors = F, sep = "\t",
                  header = F) %>%
  rename(., ID = V1) %>%
  rename(., sample_name = V2)
meta3 <- read.csv("../data/pre_label_fix/new_hiloID_link_from_Anne_8.7.18.txt", stringsAsFactors = F, sep = "\t") %>%
  rename(., ID = Library.name) %>%
  bind_rows(meta2, .) %>% # add in samples 143-160 to these samples >160
  separate(data = ., col = sample_name, sep = "_", c("popN", "family"), extra = "merge") %>%
  mutate(., n = as.integer(substr(ID, 5, 100))) %>%
  bind_rows(meta1_fixed, .) %>%
  select(., c("n", "ID", "popN", "family")) %>% # put in desired order
  mutate(., popN = as.numeric(popN)) %>%
  mutate(., zea = ifelse(popN >= 100, "maize", "mexicana")) %>%
  mutate(., symp_allo = ifelse(popN %in% c(20, 22, 33), "allopatric", "sympatric"))
# add in germplasm data from JRI for all projects ("riplasm"). 
# RIMMA is hilo maize and RIMME is hilo mex. This just gives me metadata for each population..(e.g. lat/long)
riplasm <- read.csv("../data/riplasm/riplasm.csv", stringsAsFactors = F, header = T) %>%
  mutate(., prefix = substr(RI_ACCESSION, 1, 5)) %>%
  mutate(., ind = as.integer(substr(RI_ACCESSION, 6, 100))) %>%
  dplyr::filter(., prefix %in% c("RIMMA", "RIMME")) 

# size of reference genome reads are mapped to
ref_genome_size <- sum(read.table("../data/refMaize/Zea_mays.AGPv4.dna.chr.fa.fai", 
                                  stringsAsFactors = F)$V2) # V2 is the size of each chromosome mapped to,
# including parts I won't analyze on the Pt and Mt chromosomes (but reads would pass Q filters there too)

nreads <- read.csv("../data/filtered_bam/pass1.all.metrics.Nreads", stringsAsFactors = F, 
                   header = F, sep = " ") %>%
  select(., 1:2) %>%
  rename(num_read_pass=V2) %>% # number of reads passing mapping quality filtering and deduplication
  mutate(n = as.numeric(gsub("hilo", "", V1))) %>%
  # rough coverage estimate is number of reads passing filtering * 150bp/read 
  # divided by total area they could map to in maize reference genome v4
  mutate(est_coverage = num_read_pass*150/ref_genome_size) %>%
  select(., c("n", "num_read_pass", "est_coverage"))
  
# sticking it all together
hilo <- cbind(meta3, do.call(rbind,
                             apply(meta3, 1, function(row) riplasm[riplasm$RI_ACCESSION ==  paste0(
                               ifelse(row["zea"] == "mexicana", "RIMME", "RIMMA"), 
                               sprintf("%04d", as.integer(row["popN"]))), 
                               c("RI_ACCESSION", "GEOCTY", "LOCALITY")]))) %>%
  left_join(., nreads, by = "n")

# write file with all hilo individuals
write.table(hilo, "../data/hilo_ids.txt", row.names = F, quote = F, col.names = T, sep = "\t")

# write file with only pass1 individuals -- at least some reads aligned
pass1 <- hilo[!is.na(hilo$num_read_pass), ]
pass1$pass1_ANGSD_ID <- paste0("Ind", 0:(dim(pass1)[1] - 1))
write.table(pass1, "../data/pass1_ids.txt", row.names = F, quote = F, col.names = T, sep = "\t")


# add metadata for 16 individuals from 4 lowland populations
allo4Low <- data.frame(popN=c(0,0,0,0,-10,-10,-10,-10,-20,-20,-20,-20,-30,-30,-30,-30),
                       zea = rep("maize", 16),
                       symp_allo = rep("allopatric", 16),
                       ID = paste0("4Low", c(1:4,11:14,21:24,31:34)),
                       est_coverage = rep(2, 16), # underestimate of coverage (but ok for now)
                       GEOCTY = "Mexico",
                       LOCALITY = "Lowland_4pops",
                       stringsAsFactors = F,
                       pass1_ANGSD_ID = paste0("Ind", dim(pass1)[1]:(dim(pass1)[1] + 15)))
pass1_allo4Low <- bind_rows(pass1, allo4Low)
# write files with pass1 and allo4Low individuals
write.table(pass1_allo4Low, "../data/pass1_allo4Low_ids.txt", row.names = F, quote = F, col.names = T, sep = "\t")

# aiming for more like 2x coverage for allopatric mexicana reference panel -- 
# sending new file to Anne for possible re-sequencing
filter(hilo, zea == "mexicana" & symp_allo == "allopatric") %>% 
  mutate(., est_coverage = round(est_coverage, 3)) %>% 
  .[order(.$est_coverage), ] %>%
  write.table(., "../data/resequence_allopatric_mex_for_Anne.txt",
              row.names = F, col.names = T, quote = F)
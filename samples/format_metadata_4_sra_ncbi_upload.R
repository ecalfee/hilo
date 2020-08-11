library(dplyr)
library(tidyr)
# working directory = hilo/samples
fastq_data <- read.table("HILO_fastq_metadata.txt",
                   stringsAsFactors = F, header = T, sep = "\t")
# load population metadata
load("HILO_MAIZE55_meta.RData")
# biological sample metadata
bio_cols <- c("sample_name", "organism", "isolate", "ecotype", "dev_stage",
              "geo_loc_name", "tissue", "growth_protocol", "lat_lon", 
              "isolation_source", "population")
bio <- left_join(filter(fastq_data, !duplicated(ID)),  # some IDs have multiple fastqs but only represent 1 biological sample
                 meta, 
                 by = c("ID", "popN", "family", "ind", "RI_ACCESSION")) %>%
  mutate(sample_name = ID, 
         organism = paste("Zea mays subsp.", ifelse(zea == "maize", "mays", zea)),
         isolate = paste("family", family, "seed", ind, sep = "_"),
         ecotype = RI_ACCESSION,
         dev_stage = "seedling",
         geo_loc_name = paste0("Mexico:", LOCALITY),
         tissue = "leaf",
         growth_protocol = "greenhouse",
         lat_lon = paste(ifelse(LAT > 0, paste(LAT, "N"), paste(abs(LAT, "S"))), 
                         ifelse(LON > 0, paste(LON, "E"), paste(abs(LON), "W"))),
         isolation_source = paste("Elevation:", ELEVATION, "m"),
         population = RI_ACCESSION,
         description = "Grown from seed collected as half-sib families (unknown pollen source) from populations in Mexico") %>%
  dplyr::select(., bio_cols)

write.table(bio, "bio_attributes_ncbi_sra.txt",
            quote = F, col.names = T, row.names = F, sep = "\t")

# sequencing metadata
seq_cols <- c("sample_name", "library_ID", "title", "library_strategy",
              "library_source", "library_selection", "library_layout",
              "platform", "instrument_model", "design_description",
              "filetype", "filename", "filename2")
seq <- fastq_data %>%
  mutate(sample_name = ID, 
         # library_ID already present
         title = "Whole genome sequencing of Zea mays: maize and mexicana from Mexico",
         library_strategy = "WGS",
         library_source = "GENOMIC",
         library_selection = "RANDOM",
         library_layout = "paired",
         # platform included
         instrument_model = ifelse(instrument_model == "Illumina HiSeq 6000", "Illumina NovaSeq 6000", instrument_model),
         design_description = "Seeds collected in the field were grown in the greenhouse and DNA extracted from leaf tissue. We prepared DNA libraries using a high-throughput low-volume protocol published in Rowan et al. 2019 'Nextera Low Input, Transposase Enabled protocol' (https://doi.org/10.1534/genetics.119.302406) which includes enzymatic tagmentation and sheering followed by PCR amplification and barcoding.",
         filetype = "fastq",
         filename = paste0(library_ID, "_1.fq.gz"),
         filename2 = paste0(library_ID, "_2.fq.gz"))
seq %>%
  dplyr::select(., seq_cols) %>%
  #View(.)
  write.table(., "sra_metadata_information.txt",
              quote = F, col.names = T, row.names = F, sep = "\t")


seq_long <- seq %>%
  pivot_longer(., cols = c("fastq1", "fastq2"), 
               names_to = "which_fastq", values_to = "farm_fastq") %>%
  pivot_longer(., cols = c("filename", "filename2"),
               names_to = "which_filename", values_to = "ncbi_fastq") %>%
  filter((which_filename == "filename" & which_fastq == "fastq1") |
           (which_filename == "filename2" & which_fastq == "fastq2"))

# write to files in same order for making symlinks into folder for ncbi upload
write.table(seq_long$farm_fastq, "farm_fastq_paths_4_ncbi.txt",
              quote = F, col.names = F, row.names = F, sep = "\t")
write.table(seq_long$ncbi_fastq, "ncbi_fastq_paths_4_ncbi.txt",
            quote = F, col.names = F, row.names = F, sep = "\t")

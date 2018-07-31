# identify 282 individuals from NCBI, with SRA to download, that match ts (tropial, subtropical) group
library(dplyr)
library(tidyr)
# load SRA dictionary
sra<-read.csv("../data/maize282/maize282_SRA_SRP108889",
               header = T, stringsAsFactors = F) %>%
  tidyr::separate(., "SampleName", c("SampleSet", "Inbred_Name_NCBI"), remove = F, sep = "_") %>%
  dplyr::select(., c(1:32,45,48,49)) %>%
  dplyr::mutate(Inbred_allUpperCase = toupper(Inbred_Name_NCBI))
  
# load 282 metadata from structure analysis
meta<-read.csv("../data/maize282/282_labels.csv", stringsAsFactors = F) %>%
  dplyr::rename(., Inbred_Name_282 = Inbred) %>%
  dplyr::mutate(., Inbred_allUpperCase = toupper(Inbred_Name_282))
ts<-meta[meta$Subpopulation == "ts",]
overlap_index<-which(sra$Inbred_allUpperCase %in% ts$Inbred_allUpperCase)
length(overlap_index) # some inbred's have been sequenced multiple times
length(unique(sra$Inbred_Name_NCBI[overlap_index])) # 65 of the 73 tropicals I've found on NCBI
# note that the Flint-Garcia 2005 table has 302 inbred lines (not just the 282 set??)
# discrepancies:
sum(ts$Inbred_allUpperCase %in% sra$Inbred_allUpperCase) == length(ts$Inbred) # not equal
ts$Inbred_Name_282[!(ts$Inbred_allUpperCase %in% sra$Inbred_allUpperCase)] # some tropicals not in NCBI
sra$Inbred_Name_NCBI[!(sra$Inbred_allUpperCase %in% meta$Inbred_allUpperCase)] # some NCBI not in meta data from Flint-Garcia 2005

# merge what is found as overlap
dOverlap <- dplyr::right_join(sra, ts, by = "Inbred_allUpperCase") %>%
  dplyr::select(., -Inbred_allUpperCase)
  
write.table(dOverlap, "../data/maize282/NCBI_SRA_info_for_282_ts_n=65.txt",
         col.names = T, row.names = F, quote = F, sep = "\t")

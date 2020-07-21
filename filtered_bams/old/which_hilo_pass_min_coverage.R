#!/usr/bin/env Rscript
# what is the approximate coverage for each sample,
# counting only reads that pass mapping quality > 30?
library(dplyr)
library(tidyr)
library(ggplot2)

# size of reference genome reads are mapped to
ref_genome_size <- sum(read.table("../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.fai",
                                  stringsAsFactors = F)$V2) # V2 is the size of each chromosome mapped to,
# including parts I won't analyze on the Pt and Mt and scaffolds not assigned to chromosomes (but reads would pass Q filters there too)


reads <- read.table("metrics/flagstat/multiqc/multiqc_data/multiqc_samtools_flagstat.txt",
                    sep = "\t", header = T) %>%
  dplyr::select("Sample", "total_passed") %>%
  tidyr::separate("Sample", # parse sample name into ID and lane
                  sep = " ",
                  into = c("metrics", 
                           "div1", 
                           "flagstat",
                           "div2",
                           "lane", 
                           "div3",
                           "id")) %>%
  dplyr::select("id", "lane", "total_passed") %>%
  dplyr::filter(substr(id, 1, 4) %in% c("HILO", "MAIZ")) %>% # filter out any test bams
  dplyr::mutate(group = ifelse(lane == "Maize55", "MAIZE", "HILO")) %>% # sequenced as part of Palmar Chico 55 pop (MAIZE) or this study (HILO)
  group_by(group, id) %>%
  summarise(reads = sum(total_passed),
            est_coverage = round(reads*150/ref_genome_size, 4)) # reads are 150bp
hist(reads$est_coverage)

ggplot(reads, aes(x = est_coverage, fill = group)) +
  geom_histogram() +
  facet_wrap(~group, scales = "free_x")


# what is an approximate coverage from the maize genome?
reads %>%
  group_by(group, est_coverage > 0.1) %>%
  summarise(n = n())
reads %>%
  group_by(group, est_coverage > 0.5) %>%
  summarise(n = n())
reads %>%
  group_by(group, est_coverage > 0.05) %>%
  summarise(n = n(),
            tot_coverage = sum(est_coverage))

#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(bedr)

# script to summarise diversity (pi) for homozygous ancestry windows within high introgression peaks vs. not

# get snakemake variables
pop = snakemake@params[["pop"]]
# pop = "pop366"
zea = snakemake@params[["zea"]]
# zea = "mexicana" # which homozygous ancestry are we looking at?
# Ne = 10000
# yesno = "yes"
# prefix_all = "HILO_MAIZE55"
# file_prefix = paste0("diversity/results/pi/", prefix_all, "/Ne", Ne, "_", yesno, "Boot/HOMOZYG/", zea, "/", pop)

load(snakemake@input[["meta"]])
# load(paste0("samples/", prefix_all, "_meta.RData"))
genome_file = snakemake@input[["genome"]]
# genome_file = "data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"
bed_sites = snakemake@input[["bed_sites"]]
# bed_sites = paste0("local_ancestry/results/thinnedSNPs/", prefix_all, "/whole_genome.bed")
pi_windows_file = snakemake@input[["pi_windows"]]
# pi_windows_file = paste0(file_prefix, ".pi.windows.5000.5000.pestPG")
pop_freq_file = snakemake@input[["pop_freq"]]
# pop_freq_file = paste0("local_ancestry/results/ancestry_hmm/", prefix_all, "/Ne", Ne, "_", yesno, "Boot/anc/", pop, ".anc.freq")

txt_out = snakemake@output[["txt"]] # output summary text file
# txt_out = paste0(file_prefix, ".pi.outliers_vs_not.5000.5000.txt")

# use metadata to tell if population is a maize landrace or sympatric mexicana
sympatric_zea = meta$zea[paste0("pop", meta$popN) == pop][1]

thetas_header = c("(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)", 
                  "chr", "window_center", "wattersons_theta", "pairwise_theta",
                  "tF", "tH", "tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites")

# pi_all_chr = read.table(pi_all_chr_file, stringsAsFactors = F, header = F) %>%
#   data.table::setnames(thetas_header) %>% # nSites is the number of sites within the window (individual bp) with data
#   tidyr::separate(data = .,
#                  col = `(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)`,
#                  into = c("nothing1", "indexStart", "indexStop", "firstPos_withData", "lastPos_withData", "start", "end", "nothing2")) %>%
#   dplyr::select(chr, start, end, window_center, pairwise_theta, nSites)
# 
# pi_genomewide = with(pi_all_chr, sum(pairwise_theta)/sum(nSites)) # angsd outputs pairwise_theta (tP) as a sum across sites

pi_windows = read.table(pi_windows_file, stringsAsFactors = F, header = F) %>%
  data.table::setnames(thetas_header) %>% # nSites is the number of sites within the window (individual bp) with data
  tidyr::separate(data = .,
                  col = `(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)`,
                  into = c("nothing1", "indexStart", "indexStop", "firstPos_withData", "lastPos_withData", "start", "end", "nothing2")) %>%
  dplyr::select(chr, start, end, window_center, pairwise_theta, nSites) %>%
  arrange(chr) %>% # put in chromosomal order 1-10
  mutate(chr = as.character(chr)) %>% # make chromosome a character, not integer
  mutate(start = as.numeric(start),
         end = as.numeric(end))
#with(pi_windows, sum(pairwise_theta)/sum(nSites)) #same as pi_genomewide

# combine ancestry tracts with pop ancestry frequency
# and mark outliers
freqs <- read.table(bed_sites, header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("chr", "start", "end", "pos")) %>%
  mutate(anc = read.table(pop_freq_file, header = F)$V1) %>%
  # is this tract within a high introgression outlier?
  # anc is always mexicana ancestry
  # so this is true if ancestry is high and the population is maize
  # or ancestry is low and the population is mexicana
  mutate(top_sd2 = (anc > mean(anc) + 2*sd(anc) & sympatric_zea == "maize") |
           (anc < mean(anc) - 2*sd(anc) & sympatric_zea == "mexicana")) %>%
  mutate(top_sd2 = as.numeric(top_sd2)) %>% # make TRUE = 1, FALSE = 0
  mutate(chr = as.character(chr))

#peaks = filter(freqs, top_sd2 == 1) %>%
#  bedr.merge.region(check.chr = F)

# which windows are within high introgression peaks?
# map 'top_sd2' onto windows
window_peaks = bedr(
  input = list(a = pi_windows, b = freqs), 
  method = "map", 
  check.chr = F,
  params = "-sorted -header -c 6 -o sum"
) %>%
  data.table::setnames(c(colnames(pi_windows), "n_overlap_peaks")) %>%
  dplyr::mutate(start = as.numeric(start),
                end = as.numeric(end),
                window_center = as.numeric(window_center),
                pairwise_theta = as.numeric(pairwise_theta),
                nSites = as.numeric(nSites),
                top_2sd = n_overlap_peaks > 0)


# summarise pi across windows, split by whether window overlaps peak or not (and whole genome = all windows)
window_peaks %>%
  filter(nSites > 0) %>%
  dplyr::group_by(top_2sd) %>%
  summarise(pi = sum(pairwise_theta)/sum(nSites),
            n_windows = n(),
            n_sites = sum(nSites),
            n_chr = length(unique(chr))) %>%
  mutate(windows = ifelse(top_2sd, "outlier", "non-outlier")) %>%
  bind_rows(., filter(window_peaks, nSites > 0) %>%
              summarise(pi = sum(pairwise_theta)/sum(nSites),
                        n_windows = n(),
                        n_sites = sum(nSites),
                        n_chr = length(unique(chr))) %>%
              mutate(windows = "genomewide")) %>%
  mutate(pop = pop,
         ancestry = zea) %>%
  dplyr::select(pop, ancestry, windows, pi, n_windows, n_sites, n_chr) %>%
  write.table(x = ., file = txt_out, quote = F, sep = "\t", row.names = F, col.names = T)

# option 1: for each population, plot mean pi within peaks vs. outside of peaks
# window_peaks %>%
#   filter(nSites > 0) %>%
#   mutate(pi = pairwise_theta/nSites) %>%
#   ggplot(aes(x = top_2sd, y = pi, color = top_2sd)) +
#   geom_boxplot() +
#   geom_point(data = window_peaks %>%
#                filter(nSites > 0) %>%
#                dplyr::group_by(top_2sd) %>%
#                summarise(pi = mean(pairwise_theta)/mean(nSites)))

# option 2: find peaks of a minimum length, e.g. 100kb, and
# plot pi across those peaks
# make contiguous regions out of windows within peaks:
# merged_peaks <- bedr.merge.region(x = dplyr::filter(window_peaks, top_2sd),
#                                   check.chr = F) %>%
#   mutate(length = end - start) %>%
#   arrange(as.numeric(chr))
# merged_peaks_w_pi <- bedr(
#     input = list(a = merged_peaks, b = pi_windows), 
#     method = "map", 
#     check.chr = F,
#     params = paste("-g", genome_file, "-sorted -header -c 5,6 -o sum")
#   ) %>%
#   data.table::setnames(c(colnames(merged_peaks), colnames(pi_windows)[5:6])) %>%
#   dplyr::mutate(start = as.numeric(start),
#                 end = as.numeric(end),
#                 length = as.numeric(length),
#                 pairwise_theta = as.numeric(pairwise_theta),
#                 nSites = as.numeric(nSites),
#                 pi = pairwise_theta/nSites)
#   
# # try to plot pi by length first for these contiguous outlier regions
# merged_peaks_w_pi %>%
#   filter(nSites >= 100) %>%
#   ggplot(aes(x = log(length), y = pi)) +
#   geom_point() +
#   geom_smooth()
# 
# merged_peaks_w_pi %>%
#   group_by(length) %>%
#   summarise(pi = mean(pi, na.rm = T),
#             n = n())

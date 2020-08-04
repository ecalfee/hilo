#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)

# print working directory
print(paste("R working directory:", getwd()))

# load variables from Snakefile
rmap_file = snakemake@input[["rmap"]] # recombination map
# rmap_file = "../data/linkage_map/ogut_fifthcM_map_agpv4_INCLUDE.txt"
rmap_file = "../data/linkage_map/ogut_fifthcM_map_agpv4_EXTENDED.txt"

genome_file = snakemake@input[["chr_lengths"]]
# genome_file = "../data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"
out_file = snakemake@output[["windows"]]
# out_file = "results/map_pos_1cM_windows.txt"



# get recombination map
rmap_ext <- read.table(rmap_file, stringsAsFactors = F, header = T) %>%
  data.table::setnames(c("snp", "marker_n", "cM_pos", "chr", "pos"))
rmap <- read.table("../data/linkage_map/ogut_fifthcM_map_agpv4_INCLUDE.txt", stringsAsFactors = F, header = F) %>%
  data.table::setnames(c("snp", "marker_n", "cM_pos", "chr", "pos"))



# get chromosome lengths
chr_lengths <- read.table(genome_file, stringsAsFactors = F, header = F) %>%
  data.table::setnames(c("chr", "length"))
  
# function to break a chromosome into 1cM windows
break_chr_into_1cM_windows = function(chr, rmap_chr, length){ # a chromosome, recombination map for that chr, length of that chr
  # we can interpolate/extrapolate map positions (cM) for all sites (bp) using a monotone cubic spline
  # here we want the inverse of that spline (cM -> bp), but it's not easily invertible
  # so instead we find the cM position for every bp on the chromosome
  all_cM = stats::spline(x = rmap_chr$pos, 
                         y = rmap_chr$cM_pos,
                         xmin = 1,
                         xmax = length, 
                         n = length,
                         method = "hyman")
  # then sort bp positions into windows of size 1cM
  winds = data.frame(pos = all_cM$x, pos_cM = all_cM$y,
                     stringsAsFactors = F) %>%
    mutate(window = floor(pos_cM - min(pos_cM))) %>%
    filter(!duplicated(window)) # keep only 1st position of each window
  # and put the results into bedtools format
  bed_winds = winds %>%
    mutate(start = pos - 1, start_cM = pos_cM) %>% # start of windows (bedtools is index 0)
    mutate(end = c(start[-1], tail(all_cM$x, n = 1)), # end of windows (add last bp of chromosome)
           end_cM = c(start_cM[-1], tail(all_cM$y, n = 1)),
           width = end - start,
           width_cM = end_cM - start_cM, # note: the last window will go to the end of the chr (< 1cM)
           cM_Mb = width_cM/width*10^6, # mean recomb rate in window
           chr = chr)
  return(bed_winds) 
}


test8_ext = stats::spline(x = rmap_ext$pos[rmap_ext$chr == 8], 
                       y = rmap_ext$cM_pos[rmap_ext$chr == 8],
                       xout = 1:chr_lengths$length[chr_lengths$chr == 8],
                       method = "hyman")
test8 = stats::spline(x = rmap$pos[rmap$chr == 8], 
                      y = rmap$cM_pos[rmap$chr == 8],
                      xout = 1:chr_lengths$length[chr_lengths$chr == 8],
                      method = "hyman")
test8_linear = approx(x = rmap_ext$pos[rmap_ext$chr == 8], 
                      y = rmap_ext$cM_pos[rmap_ext$chr == 8],
                      xout = 1:chr_lengths$length[chr_lengths$chr == 8],
                      method = "linear")


cor(test8_ext$y, test8$y) # not very diff except at ends of course
cor(test8_ext$y, test8_linear$y)
plot(test8_ext$y[c(T, rep(F, 10^5))], 
     test8$y[c(T, rep(F, 10^5))],
     col = "orange")
lines(test8_ext$y[c(T, rep(F, 10^5))], 
     test8_linear$y[c(T, rep(F, 10^5))],
     col = "blue")
# basically the same with linear approx
plot(test8_ext$x[c(T, rep(F, 10^5))], 
     test8_ext$y[c(T, rep(F, 10^5))],
     col = "orange", type = "l")
lines(test8_linear$x[c(T, rep(F, 10^5))], 
      test8_linear$y[c(T, rep(F, 10^5))],
      col = "blue")
# trying again:
plot(test8_ext$x[c(T, rep(F, 10^4))], 
     test8_ext$y[c(T, rep(F, 10^4))],
     col = "orange", type = "l", xlim = c(10^7, 1.2*10^7), ylim = c(20,40))
lines(test8_linear$x[c(T, rep(F, 10^4))], 
      test8_linear$y[c(T, rep(F, 10^4))],
      col = "blue")

stats::spline(x = rmap_chr$pos, 
                      y = rmap_chr$cM_pos,
                      xout = c(1, chr_lengths$length[chr_lengths$chr == 8]),
                      method = "natural")
head(rmap_chr, 1)

test8b = stats::spline(x = rmap_chr$pos, 
                      y = rmap_chr$cM_pos,
                      xout = 1:chr_lengths$length[chr_lengths$chr == 8],
                      method = "natural")
plot(test8b$x[1:100], test8b$y[1:100])
plot(test8b$x[c(T, rep(F, 10^6))], test8b$y[c(T, rep(F, 10^6))])
plot(test8$x[c(T, rep(F, 10^5))], 
     test8$y[c(T, rep(F, 10^5))],
     col = "blue", xlim = c(0, 10^7))

with(test8b[c(T, rep(F, 10^7))])

ggplot(rmap, aes(x = pos, y = cM_pos, color = as.factor(chr))) +
  geom_point() + 
  facet_wrap(~chr)

bed8 = break_chr_into_1cM_windows(chr = 8, rmap = filter(rmap, chr == 8), 
                                  length = chr_lengths$length[chr_lengths$chr == 8])
rmap_chr = filter(rmap, chr == 8)
chr = 8
length = chr_lengths$length[chr_lengths$chr == 8]
summary(bed8$cM_Mb)

plot(head(all_cM$x), head(all_cM$y))
plot(head(test8$x), head(test8$y))
# break genome into 1cM windows, chr by chr
bed = do.call(rbind, lapply(1:10, function(i)
                            break_chr_into_1cM_windows(rmap = rmap, chr = i))) %>%
  mutate(w = 1:nrow(.))

# which recombination rate quintile of the genome does a window fall into? 
# note: this is not the same as taking, e.g. the lower 5th of all windows (there are more windows per bp in high r regions of the genome)
sorted_1cM_bins <- bed %>%
  arrange(cM_Mb) %>%
  mutate(cum_pos_bp = cumsum(length)) %>%
  mutate(cum_percentile_bp = cum_pos_bp/max(cum_pos_bp))
breaks_5r_1cM <- sapply(seq(0, 1, by = .2), function(x) sorted_1cM_bins$cM_Mb[first(which(sorted_1cM_bins$cum_percentile_bp >= x))])



# .75cM/Mb is around the maize genomewide mean recomb rate

  spline = splinefun(x = rmap$pos[rmap$chr == chr], 
                y = rmap$cM_pos[rmap$chr == chr],
                method = "hyman")
bkspline = splines::backSpline(all_cM)
span_pos = spline(y = rmap$pos[rmap$chr == chr], 
                 x = rmap$cM_pos[rmap$chr == chr],
                 xout = span_cM,
                 method = "hyman")$y
# get 1cM spacing on genetic scale
tot_cM = seq(from = span_cM[1], to = span_cM[2], by = 1)
# translate to pos scale
tot_pos = spline(x = rmap$cM_pos[rmap$chr == chr],
                 y = rmap$pos[rmap$chr == chr],
                 xout = tot_cM,
                 method = "hyman")$y
sites$rpos = spline(x = rmap$pos, y = rmap$pos_cM, 
                    xout = sites$pos, method = "hyman")$y
  
  
# (2) get start and end positions for windows and each bin within windows
# (3) calculate map rate for each .2cM bin within windows
# (4) categorize map rates 1-5
# make bed files for low to high recombination: recomb_ogut_1.bed recomb_ogut_2.bed
g2 <- head(ogut)
approx(x = g2$cM_pos, y = g2$bp_pos, xout = c(-4.8, -4.0, -3.8))
# take floor to get starting position = ending position next window
# for bedtools, subtract 1 bp from starting position to get start
chr_range <- ogut %>%
  group_by(chr) %>%
  summarize(start = min(cM_pos),
            end = max(cM_pos))
window_breaks <- do.call(rbind, # calculate 1cM window breaks
        lapply(1:10, function(i){
  new_pos = approx(x = ogut$cM_pos[ogut$chr==i], # load x and y for map 
                  y = ogut$bp_pos[ogut$chr==i], 
                  xout = seq(chr_range$start[chr_range$chr==i], 
                             chr_range$end[chr_range$chr==i], 
                             by = 1))
  data.frame(chr = i,
             cM_pos = new_pos$x,
             bp_pos = floor(new_pos$y))
  })) %>%
  mutate(window_n = 1:nrow(.))
bins <- full_join(ogut, window_breaks, by = c("chr", "cM_pos")) %>%
  mutate(bp_pos = ifelse(is.na(bp_pos.x), bp_pos.y, bp_pos.x))
window_n2 <- do.call(rbind,
                     lapply(1:10, function(i)
                            as.data.frame(approx(x = window_breaks$cM_pos[window_breaks$chr == i], 
                                   y = window_breaks$window_n[window_breaks$chr == i],
                         method = "constant",
                         xout = bins$cM_pos[bins$chr == i])) %>%
                           mutate(chr = i) %>%
                           rename(cM_pos = x) %>%
                           rename(window_n = y)))
#bins2 <- 

approx(x = window_breaks$cM_pos[window_breaks$chr == 1], 
       y = window_breaks$window_n[window_breaks$chr == 1],
       method = "linear",
       xout = c(-4.8,-4,0))


# just look at head and tail of variant sites + rpos files to find non-monotonic map positions :(
vs <- bind_rows(read.table("../variant_sites/results/HILO_MAIZE55/head_regions.var.sites",
                 header = F, stringsAsFactors = F) %>%
  data.table::setnames(c("chr", "pos", "A", "a")) %>%
  mutate(rpos = read.table("../variant_sites/results/HILO_MAIZE55/head_regions.rpos",
                           header = F, stringsAsFactors = F)$V1),
  read.table("../variant_sites/results/HILO_MAIZE55/tail_regions.var.sites",
             header = F, stringsAsFactors = F) %>%
    data.table::setnames(c("chr", "pos", "A", "a")) %>%
    mutate(rpos = read.table("../variant_sites/results/HILO_MAIZE55/tail_regions.rpos",
                             header = F, stringsAsFactors = F)$V1)) %>%
  arrange(chr, pos)
table(diff(vs$rpos) < 0 & diff(vs$chr) == 0)
filter(vs, c(diff(rpos) < 0 & diff(chr) == 0, F)) %>% View()



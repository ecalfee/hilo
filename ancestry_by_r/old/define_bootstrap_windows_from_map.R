library(dplyr)
library(ggplot2)
library(bedr)

rmap = read.table("../data/linkage_map/ogut_fifthcM_map_agpv4_INCLUDE.txt", 
                  stringsAsFactors = F, header = F)
colnames(rmap) = c("SNP", "marker", "pos_cM", "chr", "pos_bp")

# round to 1 decimal place because map cM positions are
# slightly off in rounding, e.g. 0.199999
rmap$pos_cM <- round(rmap$pos_cM, 1)

# how clustered are low and high r regions?
rmap$width_bp <- c(diff(rmap$pos_bp, 1), NA)
rmap$width_cM <- round(c(diff(rmap$pos_cM, 1), NA), 1)
rmap$cM_Mb <- rmap$width_cM/(rmap$width_bp/10^6)
rmap$original <- rmap$width_cM == 0.2 # original endpoints for this .2cM segment

# 1 NA per chromosome, past final marker
rmap[c(diff(rmap$chr, 1), 1) != 0, c("width_bp", "width_cM", "cM_Mb")] <- NA
sum(is.na(rmap$cM_Mb)) == 10
sum(is.na(rmap$width_bp)) == 10
sum(is.na(rmap$width_cM)) == 10

# most markers have expected 0.2 spacing, but some larger
table(rmap$width_cM)

rmap %>%
  ggplot(aes(x = pos_cM, y = cM_Mb)) +
  geom_point() +
  facet_wrap(~chr)
rmap %>%
  filter(., chr==3) %>%
  ggplot(aes(x = pos_cM, y = cM_Mb, color = width_cM)) +
  geom_point()
rmap %>%
  ggplot(aes(x = width_cM, y = log(cM_Mb))) +
  geom_point() +
  facet_wrap(~chr)
with(rmap[complete.cases(rmap), ], cor(width_cM, cM_Mb))
with(rmap[complete.cases(rmap) & rmap$width_cM <= .4, ], cor(width_cM, cM_Mb))

# I don't think there should be higher recombination
# rates between regions with .4 vs. .2 spacing
# between markers.. my interpretation is that
# these gaps are associated with exagerrated map rates
rmap %>%
  filter(!is.na(cM_Mb)) %>%
  with(., cor(width_cM, cM_Mb))
# even when excluding the largest gaps.
# though it lessens the correlation by half
rmap %>%
  filter(!is.na(cM_Mb)) %>%
  filter(width_cM <= .4) %>%
  with(., cor(width_cM, cM_Mb))
# setting gaps to zero map distance
# is at least a slight over correction..
rmap %>%
  filter(!is.na(cM_Mb)) %>%
  filter(width_cM <= .4) %>%
  mutate(new_cM_Mb = ifelse(width_cM == .4, cM_Mb/2, cM_Mb)) %>%
  with(., cor(width_cM, new_cM_Mb))

rmap %>%
  ggplot(aes(y = log(cM_Mb), x = width_cM, color = as.factor(chr))) +
  geom_point()

# what % of the genome is 'off the map' or 'in gaps'?
# what % of SNPs are 'off the map' or 'in gaps'?
# what happens to NGSadmix estimates if I exclude these SNPs? (before bootstrapping)

# get lengths each chromosome
genome <- read.table("../data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths",
                     stringsAsFactors = F, header = F, sep = "\t") %>%
  data.table::setnames(., c("chr", "length")) %>%
  mutate(start_bp = 1)


# add intermediate 0.2 cM markers. indicate which are in 'gaps'
# or 'edges' of chromosomes
# so I can run analyses w/ and w/out 'edges'/'gaps'
genome$start_map_bp <- sapply(1:10, function(x) rmap$pos_bp[first(which(rmap$chr == x))])
genome$start_map_cM <- sapply(1:10, function(x) rmap$pos_cM[first(which(rmap$chr == x))])
genome$start_map_cM_Mb <- sapply(1:10, function(x) rmap$cM_Mb[first(which(rmap$chr == x))])
genome$end_map_bp <- sapply(1:10, function(x) rmap$pos_bp[last(which(rmap$chr == x))]) # last marker bp position
genome$end_map_cM <- sapply(1:10, function(x) rmap$pos_cM[last(which(rmap$chr == x))]) # last marker cM position
genome$end_map_cM_Mb <- sapply(1:10, function(x) rmap$cM_Mb[first(tail(which(rmap$chr == x), 2))]) # r rate between last and second to last marker
with(genome, (end_map_bp - start_map_bp)/length) # most of the chromosome is already on the map
with(genome, (start_map_cM - start_map_cM_Mb*(start_map_bp - 1)/10^6))
with(genome, (end_map_cM + end_map_cM_Mb*(length - end_map_bp)/10^6))

# write out bed file defining all the 0.2cM bins 
# and their recombination rate quintiles
# just load the extended map (goes to ends of chr's)
rmap_ext = read.table("../data/linkage_map/ogut_fifthcM_map_agpv4_EXTENDED.txt", 
                  stringsAsFactors = F, header = T)


map_pos_0.2cM_v0 <- do.call(rbind,
                          lapply(1:10, # for each chromosome
       function(x) data.frame(chr = x,
                              # get spaced cM positions each 0.2cM
                              # starting with the rounded up cM of bp=1
         pos_cM = seq(ceiling(5*rmap_ext$pos_cM[first(which(rmap_ext$chr == x))])/5,
                       rmap_ext$pos_cM[last(which(rmap_ext$chr == x))],
                       by= .2)) %>%
         # get position in bp using linear approximation
         mutate(pos_bp = round(approx(x = rmap_ext$pos_cM[rmap_ext$chr == x],
                      y = rmap_ext$pos_bp[rmap_ext$chr == x],
                      xout = pos_cM,
                      method = "linear")$y, 0)) %>%
         mutate(start = pos_bp - 1) %>% # bed coordinates start at 0
         mutate(end = c(.$start[-1], NA)) %>% # make non-overlapping by ending 1bp short of next start
         mutate(width_bp = end - start) %>%
         mutate(cM_Mb = .2/(width_bp*10^-6)) %>%
         filter(!is.na(end)) # last position just marks end of the last interval, not it's own start
       )) %>% 
  #name intervals, e.g. w9
  mutate(window = paste0("w", 1:nrow(.))) %>%
  left_join(., rmap[ , c("chr", "pos_cM", "original")], by = c("chr", "pos_cM")) %>%
  # mark if interval was one of the original 0.2cM windows mapped w/ 2 original markers
  mutate(original_map0.2markers = sapply(original, isTRUE))

# also get 1cM windows and their mean recombination rates
map_pos_1cM_v0 <- do.call(rbind,
                       lapply(1:10, # for each chromosome
                              function(x) data.frame(chr = x,
                                                     # get spaced cM positions each 0.2cM
                                                     # starting with the rounded up cM of bp=1
                                                     pos_cM = seq(ceiling(5*rmap_ext$pos_cM[first(which(rmap_ext$chr == x))])/5,
                                                                  rmap_ext$pos_cM[last(which(rmap_ext$chr == x))],
                                                                  by = 1)) %>%
                                # get position in bp using linear approximation
                                mutate(pos_bp = round(approx(x = rmap_ext$pos_cM[rmap_ext$chr == x],
                                                             y = rmap_ext$pos_bp[rmap_ext$chr == x],
                                                             xout = pos_cM,
                                                             method = "linear")$y, 0)) %>%
                                mutate(start = pos_bp - 1) %>% # bed coordinates start at 0
                                mutate(end = c(.$start[-1], NA)) %>% # make non-overlapping by ending 1bp short of next start
                                mutate(width_bp = end - start) %>%
                                mutate(cM_Mb = 1/(width_bp*10^-6)) %>%
                                filter(!is.na(end)) # last position just marks end of the last interval, not it's own start
                       )) %>% 
  #name intervals, e.g. w9
  mutate(window = paste0("W", 1:nrow(.))) 


# ooops, top 5th of cM windows is a lot smaller than top 5th of genome
# I want to define quantiles from physical space covered, not proportion windows

# get quantile bin for recombination rate 
sorted_0.2cM_bins <- map_pos_0.2cM_v0 %>%
    arrange(cM_Mb) %>%
    mutate(cum_pos_bp = cumsum(width_bp)) %>%
    mutate(cum_percentile_bp = cum_pos_bp/max(cum_pos_bp))
sorted_1cM_bins <- map_pos_1cM_v0 %>%
  arrange(cM_Mb) %>%
  mutate(cum_pos_bp = cumsum(width_bp)) %>%
  mutate(cum_percentile_bp = cum_pos_bp/max(cum_pos_bp))
breaks_5r_0.2cM <- sapply(seq(0, 1, by = .2), function(x) sorted_0.2cM_bins$cM_Mb[first(which(sorted_0.2cM_bins$cum_percentile_bp >= x))])
breaks_5r_1cM <- sapply(seq(0, 1, by = .2), function(x) sorted_1cM_bins$cM_Mb[first(which(sorted_1cM_bins$cum_percentile_bp >= x))])
  
# label bins in data by their quantile  
map_pos_0.2cM <- map_pos_0.2cM_v0 %>%
  mutate(bin_r5 = cut(cM_Mb, 
                      breaks = breaks_5r_0.2cM, 
                      right = T,
                      include.lowest = T)) %>%
  dplyr::select(., c("chr", "start", "end", "window", "cM_Mb", "bin_r5", "pos_cM", "original_map0.2markers"))
  
map_pos_1cM <- map_pos_1cM_v0 %>%
  mutate(bin_r5 = cut(cM_Mb, 
                      breaks = breaks_5r_1cM, 
                      right = T,
                      include.lowest = T)) %>%
  dplyr::select(., c("chr", "start", "end", "window", "cM_Mb", "bin_r5", "pos_cM"))

# a little less than half the intervals have 2 original mapped markers:
table(map_pos_0.2cM$original_map0.2markers)
table(map_pos_0.2cM$bin_r5)
table(map_pos_1cM$bin_r5) # similar story but fewer extremes
# hist quantile differences neighboring windows
hist(diff(as.integer(map_pos_0.2cM$bin_r5)))
hist(diff(as.integer(map_pos_1cM$bin_r5)))
# how spaced are they? can I merge into 1cM regions? 
# worth doing this broader scale in addition -- it's mostly a centromere/telomere thing
map_pos_0.2cM %>%
  ggplot(., aes(x = pos_cM, y = log(cM_Mb), color = original_map0.2markers)) +
  geom_point() +
  ggtitle("0.2cM recombination windows") +
  facet_wrap(~chr)
ggsave("plots/map_pos_0.2cM_windows_log_recomb_rate.png",
       width = 8, height = 8, units = "in")
map_pos_1cM %>%
  ggplot(., aes(x = pos_cM, y = log(cM_Mb))) +
  geom_point() +
  ggtitle("1cM recombination windows") +
  facet_wrap(~chr)
ggsave("plots/map_pos_1cM_windows_log_recomb_rate.png",
       width = 8, height = 8, units = "in")

# print lists for each quantile
write.table(map_pos_0.2cM, "results/map_pos_0.2cM_windows.txt",
            quote = F, col.names = T, row.names = F, sep = "\t")
write.table(map_pos_1cM, "results/map_pos_1cM_windows.txt",
            quote = F, col.names = T, row.names = F, sep = "\t")

# I may need to exclude windows with < 5 SNPs or something similar. First get SNPs in windows.
r_windows_excl <- c("W382") # for now just exclude W382 for having no SNPs

# make 100 bootstrap replicate lists of windows:
set.seed(100) # set seed for replication
for (r in 1:5){ # 5 recombination rate quintiles
  r_windows <- map_pos_1cM$window[map_pos_1cM$bin_r5 == levels(map_pos_1cM$bin_r5)[r]]
  r_windows_incl <- r_windows[!(r_windows %in% r_windows_excl)]
  # boot0 is the original data/windows:
  write.table(x = r_windows_incl,
              file = paste0("results/bootstrap/windows_1cM/r5_recomb", r, "/boot", 0, ".list"),
              quote = F, col.names = F, row.names = F, sep = "\t")
  for (b in 1:100){ # 100 bootstrap samples 
    # (same size = same # windows as original data 
    # but windows are sampled w/ replacement)
    write.table(x = sample(x = r_windows_incl,
                         size = length(r_windows_incl),
                         replace = T),
                file = paste0("results/bootstrap/windows_1cM/r5_recomb", r, "/boot", b, ".list"),
                quote = F, col.names = F, row.names = F, sep = "\t")
  }
}


# extra stuff:
windows_10kb <- read.table("../data/refMaize/windows_10kb/whole_genome.bed",
                           header = F, sep = "\t", stringsAsFactors = F) %>%
  data.table::setnames(c("chr", "start", "end", "N_window_10kb"))%>%
  mutate(chr = as.character(chr))
windows_1kb <- read.table("../data/refMaize/windows_1kb_whole_genome.bed",
                           header = F, sep = "\t", stringsAsFactors = F) %>%
  data.table::setnames(c("chr", "start", "end"))%>%
  mutate(chr = as.character(chr))
r_0.2cM_by_10kb_window <- bedr(
  engine = "bedtools", 
  input = list(a = windows_10kb, b = mutate(map_pos_0.2cM, chr = as.character(chr))), 
  method = "map", 
  params = "-c 5 -o mean -g ../data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths",
  check.chr = F
) %>%
  data.table::setnames(c("chr", "start", "end", "N_window_10kb", "mean_cM_Mb")) %>%
  mutate(mean_cM_Mb = as.numeric(mean_cM_Mb))
r_0.2cM_by_1kb_window <- bedr(
  engine = "bedtools", 
  input = list(a = windows_1kb, b = mutate(map_pos_0.2cM, chr = as.character(chr))), 
  method = "map", 
  params = "-c 5 -o mean -g ../data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths",
  check.chr = F
) %>%
  data.table::setnames(c("chr", "start", "end", "mean_cM_Mb")) %>%
  mutate(mean_cM_Mb = as.numeric(mean_cM_Mb))

# quantiles over 10kb windows pretty much match quantiles over bp
# because these windows are very short compared to changes in recomb. rate
# but it's more accurate to get quantiles from ind. bp as I do above
# using cumsum() function (this is best)
round(quantile(r_0.2cM_by_10kb_window$mean_cM_Mb, na.rm = T), 2)
round(quantile(r_0.2cM_by_1kb_window$mean_cM_Mb, na.rm = T), 2)


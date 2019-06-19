ogut <- read.table("../data/linkage_map/ogut_fifthcM_map_agpv4_INCLUDE.txt",
                   stringsAsFactors = F, header = F)
colnames(ogut) <- c("snp", "marker_n", "cM_pos", "chr", "bp_pos")
# (1) break genome into 1cM windows, starting with lowest marker per chromosome
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

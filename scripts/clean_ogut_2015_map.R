#!/usr/bin/env Rscript
# this script drops 'bad' SNPs from the 0.2 cM Ogut 2015 maize 
# recombination map where a few such SNPs appear to be out
# of order, due to map or assembly error, e.g. @ 8.8 on chr2 below:
#S_3403255    M1114    8.4    2    3432355
#S_3432232    M1115    8.6    2    3461331
#S_3461208    M1116    8.8    2    3284078
#S_3490185    M1117    9    2    3492426
#S_3519161    M1118    9.2    2    3521396
# notably, this affects a few larger contiguos regions on chr 7
library(dplyr)
library(ggplot2)


rmap = read.table("../data/linkage_map/ogut_fifthcM_map_agpv4.txt", 
                stringsAsFactors = F, header = F)
colnames(rmap) = c("SNP", "marker", "pos_cM", "chr", "pos_bp")
rmap = rmap[!is.na(rmap$SNP),]

# find positions that have a chromosome that matches neither the chrom in front nor behind them
i_missChr = which(lead(rmap$chr) != rmap$chr & lag(rmap$chr) != rmap$chr)
# remove these SNPs
rmap_1 = rmap[-i_missChr,]
# confirmed no more
length(which(lead(rmap_1$chr) != rmap_1$chr & lag(rmap_1$chr) != rmap_1$chr)) == 0
# difference in chromosome number between current and next mapped position
table(lead(rmap_1$chr) - rmap_1$chr)
i_missChr2 = which(!(lead(rmap_1$chr) - rmap_1$chr) %in% 0:1)
# clearly a small segment of chr7 stuck within chr2 -- I'll remove by hand
#lapply(i_missChr2, function(i) rmap_1[(i-3):(i+3), ])
rmap_1 = filter(rmap_1, !(marker %in% c("M1786", "M1787", "M1788")))
table(lead(rmap_1$chr) - rmap_1$chr) # looks good! All markers out of order due to chromosome are now removed

# now I look for SNPs embedded within chromosomes that appear to be out-of-order
# afterwards I deal with the case of SNPs near the chromosome transition boundary

# find SNPs that are on the same chromosome and greater than all the SNPs up to 
#3 behind and 3 ahead of themselves on the map
tooBig3 = function(map) {
  map$chr == lag(map$chr) &
    map$chr == lag(lag(map$chr)) &
    map$chr == lag(lag(lag(map$chr))) &
    map$chr == lead(map$chr) &
    map$chr == lead(lead(map$chr)) &
    map$chr == lead(lead(lead(map$chr))) &
    map$pos_bp > lag(map$pos_bp) & 
    map$pos_bp > lag(lag(map$pos_bp)) & 
    map$pos_bp > lag(lag(lag(map$pos_bp))) &
    map$pos_bp > lead(map$pos_bp) & 
    map$pos_bp > lead(lead(map$pos_bp)) & 
    map$pos_bp > lead(lead(lead(map$pos_bp)))
}

# find SNPs that are on the same chromosome and are less than all the SNPs up to 
# 3 behind and 3 ahead of themselves on the map
tooSmall3 = function(map) {
  map$chr == lag(map$chr) &
    map$chr == lag(lag(map$chr)) &
    map$chr == lag(lag(lag(map$chr))) &
    map$chr == lead(map$chr) &
    map$chr == lead(lead(map$chr)) &
    map$chr == lead(lead(lead(map$chr))) &
  map$pos_bp < lag(map$pos_bp) & 
    map$pos_bp < lag(lag(map$pos_bp)) & 
    map$pos_bp < lag(lag(lag(map$pos_bp))) &
    map$pos_bp < lead(map$pos_bp) & 
    map$pos_bp < lead(lead(map$pos_bp)) & 
    map$pos_bp < lead(lead(lead(map$pos_bp)))
}

# find SNPs that are on the same chromosome and greater than all the SNPs up to 
#2 behind and 2 ahead of themselves on the map
tooBig2 = function(map) {
  map$chr == lag(map$chr) &
    map$chr == lag(lag(map$chr)) &
    map$chr == lead(map$chr) &
    map$chr == lead(lead(map$chr)) &
    map$pos_bp > lag(map$pos_bp) & 
    map$pos_bp > lag(lag(map$pos_bp)) & 
    map$pos_bp > lead(map$pos_bp) & 
    map$pos_bp > lead(lead(map$pos_bp))
}

# find SNPs that are on the same chromosome and are less than all the SNPs up to 
# 2 behind and 2 ahead of themselves on the map
tooSmall2 = function(map) {
  map$chr == lag(map$chr) &
    map$chr == lag(lag(map$chr)) &
    map$chr == lead(map$chr) &
    map$chr == lead(lead(map$chr)) &
    map$pos_bp < lag(map$pos_bp) & 
    map$pos_bp < lag(lag(map$pos_bp)) & 
    map$pos_bp < lead(map$pos_bp) & 
    map$pos_bp < lead(lead(map$pos_bp))
}

# find SNPs that are on the same chromosome and greater than all the SNPs up to 
#1 behind and 1 ahead of themselves on the map
tooBig1 = function(map) {
  map$chr == lag(map$chr) &
    map$chr == lead(map$chr) &
    map$pos_bp > lag(map$pos_bp) & 
    map$pos_bp > lead(map$pos_bp)
}

# find SNPs that are on the same chromosome and are less than all the SNPs up to 
# 1 behind and 1 ahead of themselves on the map
tooSmall1 = function(map) {
  map$chr == lag(map$chr) &
    map$chr == lead(map$chr) &
    map$pos_bp < lag(map$pos_bp) & 
    map$pos_bp < lead(map$pos_bp)
}


# apply tooBig and tooSmall filter recursively
# until no more SNPs are filtered out
recursive_filter3 = function(map){
  nSNP = nrow(map)
  # do one step of filtering
  map2 = map %>%
  filter(., !tooBig3(.) | is.na(tooBig3(.))) %>%
    filter(., !tooSmall3(.)| is.na(tooSmall3(.)))
  if (nrow(map2) < nrow(map)){ # SNPs still being filtered out
    return(recursive_filter3(map2))
  } else{ # otherwise done filtering, return current map
    return(map2)
  }
}
recursive_filter2 = function(map){
  nSNP = nrow(map)
  # do one step of filtering
  map2 = map %>%
    filter(., !tooBig2(.) | is.na(tooBig2(.))) %>%
    filter(., !tooSmall2(.)| is.na(tooSmall2(.)))
  if (nrow(map2) < nrow(map)){ # SNPs still being filtered out
    return(recursive_filter2(map2))
  } else{ # otherwise done filtering, return current map
    return(map2)
  }
}
recursive_filter1 = function(map){
  nSNP = nrow(map)
  # do one step of filtering
  map2 = map %>%
    filter(., !tooBig1(.) | is.na(tooBig1(.))) %>%
    filter(., !tooSmall1(.)| is.na(tooSmall1(.)))
  if (nrow(map2) < nrow(map)){ # SNPs still being filtered out
    return(recursive_filter1(map2))
  } else{ # otherwise done filtering, return current map
    return(map2)
  }
}

recursive_filter_edge = function(map){
  nSNP = nrow(map)
  # do one step of filtering to find any remaining SNPs that are smaller than the next SNP
  # but still on the same chromosome
  bad_edge = which(lead(map$pos_bp) <= map$pos_bp & lead(map$chr) == map$chr)
  if (length(bad_edge) > 0){ # SNPs still being filtered out
    print(paste("some edge problems:", map[bad_edge,]))
    return(recursive_filter_edge(map[-bad_edge, ])) 
    # drop SNP positions that are greater than the next's SNPs pos
  } else{ # otherwise done filtering, return current map
    return(map)
  }
}


recursive_filter_all = function(map){
  # first filter out SNPs bigger/smaller than their 3 neighbors before/after
  # then 2 neighbors
  # then 1 neighbor
  map2 = recursive_filter1(recursive_filter2(recursive_filter3(map)))
  #map2 = recursive_filter_edge(recursive_filter1(recursive_filter2(recursive_filter3(map))))
  # finally filter out any remaining SNPs that are on chromosome boundaries and out of order
} # I am not sure why the last filter for edge problems isn't working -- will need to check
  
rmap_2 = recursive_filter_all(rmap_1)
#dim(rmap_2)

# this filtering doesn't change the general recombination map/pattern
make_plots = function(rmap_1, rmap_2){

  par(mfrow=c(2,1))
  plot(as.numeric(row.names(rmap_1)), rmap_1$pos_bp, col = "blue", main = paste0("pre-filtering n=", nrow(rmap_1)))
  plot(as.numeric(row.names(rmap_2)), rmap_2$pos_bp, main = paste0("post-filtering n=", nrow(rmap_2)))
  par(mfrow=c(1,1))

  # zoom in on a pos_cM to pos_bp map. from the general shape of these maps it makes more sense
  # to extend the tail recombination rate to markers beyond the map boundaries (on end of chroms)
  # than it does to use a chromosome-average for these markers because rates are faster on the ends
  # than near the centromere, so this should be a better proxy
  rmap_2 %>%
    ggplot(., aes(pos_bp, pos_cM)) +
    geom_point() +
    facet_wrap(.~chr)
  rmap_1 %>%
    #filter(., chr == 1) %>%
    ggplot(., aes(pos_bp, pos_cM)) +
    geom_point(color = "blue") +
    facet_wrap(.~chr)  
}
# make an output file for included markers
write.table(rmap_2,
            "../data/linkage_map/ogut_fifthcM_map_agpv4_INCLUDE.txt",
            sep = "\t", col.names = F, row.names = F, quote = F)

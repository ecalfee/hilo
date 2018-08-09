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
table(lead(rmap_1$chr) - rmap_1$chr)
i_missChr2 = which(!(lead(rmap_1$chr) - rmap_1$chr) %in% 0:1)
# clearly a small segment of chr7 stuck within chr2 -- I'll remove by hand
lapply(i_missChr2, function(i) rmap_1[(i-3):(i+3), ])
rmap_1 = filter(rmap_1, !(marker %in% c("M1786", "M1787", "M1788")))
table(lead(rmap_1$chr) - rmap_1$chr) # looks good! All markers out of order due to chromosome are now removed


# difference between current bp position and next position
rmap_1$diff = lead(rmap_1$pos_bp) - rmap_1$pos_bp
# different chromosome, don't measure difference across chromosomes
#rmap_1$diff[rmap_1$chr != lag(rmap_1$chr)] <- NA


# find SNPs that are greater than 3 behind and 3 ahead of themselves
tooBig = function(map) {
  (map$pos_bp > lag(map$pos_bp) | map$chr != lag(map$chr)) & 
                   (map$pos_bp > lag(lag(map$pos_bp)) | map$chr != lag(lag(map$chr))) & 
                   (map$pos_bp > lag(lag(lag(map$pos_bp))) | map$chr != lag(lag(lag(map$chr)))) &
                   map$pos_bp > lead(map$pos_bp) & 
                   map$pos_bp > lead(lead(map$pos_bp)) & 
                   map$pos_bp > lead(lead(lead(map$pos_bp)))
}


# find SNPs that are less than 3 behind and 3 ahead of themselves
# WHY DOES THIS KEEP EATING AWAY AT MY STUFF??
tooSmall = function(map) {
  map$pos_bp < lag(map$pos_bp) & 
                   map$pos_bp < lag(lag(map$pos_bp)) & 
                   map$pos_bp < lag(lag(lag(map$pos_bp))) &
                   (map$pos_bp < lead(map$pos_bp) & map$chr == lead(map$chr)) & 
                   (map$pos_bp < lead(lead(map$pos_bp)) | map$chr == lead(lead(map$chr))) & 
                   (map$pos_bp < lead(lead(lead(map$pos_bp))) | map$chr == lead(lead(lead(map$chr))))
}

tooBig = function(map) {
  map$pos_bp > lag(map$pos_bp) & 
    map$pos_bp > lag(lag(map$pos_bp)) & 
    map$pos_bp > lag(lag(lag(map$pos_bp))) &
    map$pos_bp > lead(map$pos_bp) & 
    map$pos_bp > lead(lead(map$pos_bp)) & 
    map$pos_bp > lead(lead(lead(map$pos_bp)))
}


# find SNPs that are less than 3 behind and 3 ahead of themselves
tooSmall = function(map) {
  map$pos_bp < lag(map$pos_bp) & 
    map$pos_bp < lag(lag(map$pos_bp)) & 
    map$pos_bp < lag(lag(lag(map$pos_bp))) &
    map$pos_bp < lead(map$pos_bp) & 
    map$pos_bp < lead(lead(map$pos_bp)) & 
    map$pos_bp < lead(lead(lead(map$pos_bp)))
}


rmap_2 = rmap_1 %>%
  filter(., !tooBig(.)) %>%
  filter(., !tooSmall(.)) %>%
  filter(., !tooBig(.)) %>%
  filter(., !tooSmall(.)) %>%
  filter(., !tooBig(.)) %>%
  filter(., !tooSmall(.)) %>%
  filter(., !tooBig(.)) %>%
  filter(., !tooSmall(.)) %>%
  filter(., !tooBig(.)) %>%
  filter(., !tooSmall(.)) %>%
  filter(., !tooBig(.)) %>%
  filter(., !tooSmall(.)) %>%
  filter(., !tooBig(.)) %>%
  filter(., !tooSmall(.)) %>%
  filter(., !tooBig(.)) %>%
  filter(., !tooSmall(.)) %>%
  filter(., !tooBig(.)) %>%
  filter(., !tooSmall(.)) %>%
  filter(., !tooBig(.)) %>%
  filter(., !tooSmall(.)) %>%
  filter(., !tooBig(.)) %>%
  filter(., !tooSmall(.)) %>%
  filter(., !tooBig(.)) %>%
  filter(., !tooSmall(.)) %>%
  filter(., !tooBig(.)) %>%
  filter(., tooSmall(.))


rmap_2a = rmap_1[!i_tooBig(rmap_1),]
rmap_2b = rmap_2a[!i_tooSmall(rmap_2a),]


rmap$chr_prev = c(0, rmap$chr[1:(nrow(rmap) - 1)]) # previous chromosome
# exclude markers on chromosomes out of order
excl_candidate = which(rmap$chr < rmap$chr_prev)
# visualize each one in context to identify correctly
rmap[sort(c(excl_candidate -4, excl_candidate -3, excl_candidate - 2, excl_candidate - 1, excl_candidate, 
       excl_candidate + 1, excl_candidate + 2, excl_candidate + 3)), ]
excl_i1 = c(1624, 1640, 1641, 1642, 2631, 4114, 4169, 4508, 4589, 6085, 6245)
exclude1 = rmap[excl_i1, c("SNP", "marker", "pos_cM", "chr", "pos_bp")] # exclude
rmap1 = rmap[-excl_i1, c("SNP", "marker", "pos_cM", "chr", "pos_bp")] # keep
rmap1$chr_prev = c(0, rmap1$chr[1:(nrow(rmap1) - 1)]) # reset previous chromosome
rmap1[rmap1$chr < rmap1$chr_prev,] # none - good



#### THE NEXT STEP IS INCOMPLETE - SOME DIFF. CASES OF MULTIPLE BAD SNPs IN A ROW ETC##
# now look for markers out of order within chromosome
diff_bp = diff(rmap1$pos_bp)
rmap1$diff = c(0, diff_bp) # look 1 position back
diff_bp2 = diff(rmap1$pos_bp, lag = 2) # look 2 back
rmap1$diff2 = c(0, 0, diff_bp2)
# first eliminate positions less than the 1 previous and 2 previous bp
rmap1$filter_out = rmap1$diff < 0 & rmap1$diff2 < 0 & rmap1$chr == rmap1$chr_prev
# visualize
#rmap1[sort(c(which(rmap1$filter_out==T)-2, 
#             which(rmap1$filter_out==T)-1, 
#             which(rmap1$filter_out==T),
#             which(rmap1$filter_out==T) + 1)), ]
exclude2 = rmap1[rmap1$filter_out, c("SNP", "marker", "pos_cM", "chr", "pos_bp")]
rmap2 = rmap1[!rmap1$filter_out, c("SNP", "marker", "pos_cM", "chr", "pos_bp")]
# then eliminate positions less than the 1 previous bp
rmap2$chr_prev = c(0, rmap2$chr[1:(nrow(rmap2) - 1)]) # reset previous chromosome
diff_bp = diff(rmap2$pos_bp)
rmap2$diff = c(0, diff_bp)
exclude3 = rmap2[rmap2$diff < 0 & rmap2$chr_prev == rmap2$chr, c("SNP", "marker", "pos_cM", "chr", "pos_bp")]
rmap3 = rmap2[!(rmap2$diff < 0 & rmap2$chr_prev == rmap2$chr), c("SNP", "marker", "pos_cM", "chr", "pos_bp")]
diff_bp3 = diff(rmap3$pos_bp)
rmap3$diff = c(0, diff_bp3)
rmap3[rmap3$diff <0 & rmap3$chr_prev == rmap3$chr,]
# make an output file for included markers
write.table(rmap1[!rmap1$filter_out, 
                 c("SNP", "marker", "pos_cM", "chr", "pos_bp")],
                 "../data/linkage_map/ogut_fifthcM_map_agpv4_INCLUDE.txt",
            sep = "\t", col.names = F, row.names = F, quote = F)
# and make a separate output file for excluded markers
write.table(rbind(exclude1, rmap1[rmap1$filter_out, 
                 c("SNP", "marker", "pos_cM", "chr", "pos_bp")]),
            "../data/linkage_map/ogut_fifthcM_map_agpv4_EXCLUDE.txt",
            sep = "\t", col.names = F, row.names = F, quote = F)


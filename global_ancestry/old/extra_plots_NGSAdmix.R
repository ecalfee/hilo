# p_allo_elev <- d %>%
#   dplyr::mutate(., mexicana = anc2,
#                 est_coverage = ifelse(is.na(est_coverage), 30, est_coverage)) %>%
#   filter(., symp_allo == "allopatric") %>%
#   ggplot(., aes(x = LOCALITY, y = mexicana, color = LOCALITY,
#                 size = log(est_coverage), shape = zea)) +
#   geom_jitter() +
#   ylab("Proportion mexicana ancestry") +
#   xlab("Location") +
#   ggtitle("Clines in mexicana ancestry across elevation") +
#   theme(legend.position="bottom") +
#   labs(color = "Location", shape = "Subspecies", size = "log(Coverage)") +
#   theme_classic()
# ggsave(filename = "../test/ngsadmix_allo.png",
#        plot = p_allo_elev,
#        device = "png", height = 8, width = 12, units = "in", dpi = 300)
# 


# # plot 'STRUCTURE-like' ancestry plots
# png(paste0("../plots/NGSadmix_K", K, ".png"), # saves plot as pin in ../plots/
#     height = 5, width = 8, units = "in", res = 150)
# d %>%
#   dplyr::select(., colnames(admix)) %>%
#   t(.) %>%
#   barplot(height = .,
#           col = colorsK[1:K],
#           space=0,
#           border=NA,
#           main = "all HILO populations",
#           xlab="Individuals",
#           ylab="admixture") 
# mtext(paste("est. between-ancestry Fst:", round(Fst,2)))
# dev.off()
# 
# # now plot again with individuals with very low coverage <0.05x filtered out      
# png(paste0("../plots/NGSadmix_K", K, "_over_", min_coverage, "x_coverage.png"),
#     height = 5, width = 8, units = "in", res = 150)
# d %>%
#   filter(., est_coverage >= min_coverage) %>%
#   select(., colnames(admix)) %>%
#   t(.) %>%
#   barplot(height = .,
#           col = colorsK[1:K],
#           space=0,
#           border=NA,
#           main = "all HILO populations",
#           xlab=paste0("Individuals with est. coverage >= ", min_coverage, "x"),
#           ylab="admixture") 
# mtext(paste("est. between-ancestry Fst:", round(Fst,2)))
# dev.off()
# 
# for (g in unique(d$group)){
#   png(paste0("../plots/NGSadmix_K", K, "_", g, ".png"), # saves plot as pin in ../plots/
#       height = 5, width = 8, units = "in", res = 150)
#   d %>%
#     filter(., group == g) %>%
#     filter(., est_coverage >= min_coverage) %>%
#     select(., colnames(admix)) %>%
#     t(.) %>%
#     barplot(height = .,
#             col = colorsK[1:K],
#             space=0,
#             border=NA,
#             main = g,
#             xlab=paste0("Individuals with est. coverage >= ", min_coverage, "x"),
#             ylab="admixture")  
#   dev.off()
# }
# 
# 
# # by location
# if (K == 2){
#   p = d %>%
#     filter(., est_coverage >= .05) %>%
#     mutate(group = factor(paste(symp_allo, zea, sep = "_"), ordered = T,
#                           levels = labels_maize2mex)) %>%
#     ggplot(., aes(x = group, y = anc1, color = group)) +
#     geom_violin() + 
#     geom_point(size=1, position = position_jitter(w=0.05)) +
#     scale_colour_manual(values = colors_maize2mex) +
#     #scale_colour_manual(values = c("yellow", "darkblue", "orange", "blue")) +
#     theme_minimal() + 
#     labs(y = "proportion 'mexicana-like' ancestry")
#   p + facet_wrap(~LOCALITY) # show plot
#   # save plot
#   ggsave("plots/NGSadmix_proportion_mexicana-like_by_location.png", plot = p + facet_wrap(~LOCALITY), device = png(), 
#          width = 12, height = 8, units = "in",
#          dpi = 200)
# }
# 
# 
# # make NGSadmix plot for poster:
# # get elevation data
# meta.pops = read.table("../data/riplasm/gps_and_elevation_for_sample_sites.txt",
#                  stringsAsFactors = F, header = T, sep = "\t")
# meta.ind = left_join(d[ , c("ID", "popN", "total_reads_pass", "est_coverage", "group")], meta.pops, by = "popN") %>%
#   filter(., est_coverage >= min_coverage) %>% # only include individuals with at least .05x coverage
#   .[with(., order(group, ELEVATION)), ]
# # combined with admixture results
# d_meta.ind = left_join(dplyr::select(d, ID, starts_with("anc")), meta.ind, by = "ID")
# 
# 
# 
# 
# # now plot again with individuals with very low coverage <0.05x filtered out      
# png(paste0("plots/", PREFIX, "NGSadmix_K", K, "_over_", min_coverage, "x_coverage_wElevation_for_poster.png"),
#     height = 5, width = 8, units = "in", res = 300)
# par(mar=c(4.1,4.1,4.1,4.1))
# bar = d_meta.ind %>%
#   #dplyr::select(., colnames(admix)) %>%
#   dplyr::select(., rev(colnames(admix))) %>% # reverse order so blue = mexicana
#   t(.) %>%
#   barplot(height = .,
#           col = colorsK[1:K],
#           space=0,
#           border=NA, xaxt = "n")
# title(main = "Admixture in highland maize and mexicana",
#       ylab=paste0("Ancestry proportion K=", K, " (NGSAdmix)"))
# #title(xlab="Allopatric Maize | Allopatric Mexicana | Sympatric Maize | Sympatric Mexicana", 
# #      line = 1)
# text(x = tapply(bar, d_meta.ind$group, mean), par("usr")[3], srt = 60, adj= 1, xpd = TRUE,
#      labels = unique(d_meta.ind$group), cex=0.45)
# #title(xlab=paste("Allopatric Ref.", "Sympatric Maize", "Sympatric Mexicana",
# #      sep = "               |               "),
# #      line = 1)
# par(new = TRUE)
# plot(x = bar, y = d_meta.ind$ELEVATION, ylim = c(1000, 3000), axes = F, type = "p", cex = .1, 
#      xlab = "", ylab = "")
# axis(side = 4, at = c(1000, 2000, 3000), 
#      tick = T, labels = T, col = "black")
# mtext("Elevation (m)", side=4, line=2.5)
# par(mar=c(5.1,4.1,4.1,2.1)) # set back default
# dev.off()
# 
# # make simple plots for SNPs pruned every 1000th SNP
# file_prefix_1000 = paste0("../data/geno_lik/merged_pass1_all_alloMaize4Low_16/prunedBy1000/NGSAdmix/K", K)
# 
# admix_1000 <- read.table(paste0(file_prefix_1000, ".qopt"))
# colnames(admix_1000) <- paste0("anc", 1:K) #c("anc1", "anc2")
# admix_1000 <- admix_1000[, K:1] # reverse order for pruned by 1000
# d_1000 <- bind_cols(pass1_allo4Low, admix_1000)  %>%
#   arrange(., popN) %>%
#   arrange(., zea) %>%
#   arrange(., symp_allo) %>%
#   mutate(., group = paste(symp_allo, zea, sep = "_"))
# d_meta.ind_1000 = left_join(d_1000, meta, by = c("popN", "zea", "symp_allo", 
#                                    "RI_ACCESSION", "GEOCTY", "LOCALITY")) %>%
#   filter(., est_coverage >= .05) %>% # only include individuals with at least .05x coverage
#   .[with(., order(group, ELEVATION)), ]
# 
# 
# # plot results for SNPs pruned down to 1 per 1000 SNPs
# png(paste0("../plots/NGSadmix_K", K, "_over_0.05x_coverage_wElevation_for_poster_pruned1000.png"),
#     height = 5, width = 8, units = "in", res = 300)
# par(mar=c(4.1,4.1,4.1,4.1))
# bar_1000 = d_meta.ind_1000 %>%
#   select(., colnames(admix_1000)) %>%
#   t(.) %>%
#   barplot(height = .,
#           col = colorsK[1:K],
#           space=0,
#           border=NA, xaxt = "n")
# title(main = "Admixture in highland maize and mexicana; every 1000th SNP",
#       ylab=paste0("Ancestry proportion K=", K, " (NGSAdmix)"))
# title(xlab="Allopatric Maize | Allopatric Mexicana | Sympatric Maize | Sympatric Mexicana", 
#       line = 1)
# 
# par(new = TRUE)
# plot(x = bar_1000, y = d_meta.ind$ELEVATION, ylim = c(1000, 3000), axes = F, type = "p", cex = .1, 
#      xlab = "", ylab = "")
# axis(side = 4, at = c(1000, 2000, 3000), 
#      tick = T, labels = T, col = "black")
# mtext("Elevation (m)", side=4, line=2.5)
# par(mar=c(5.1,4.1,4.1,2.1)) # set back default
# dev.off()
# 
# # do global ancestry estimates
# # from NGSadmix confirm the ancestry ~ r pattern?
# # get ngsadmix k=2 estimates for 5 bins of recombination rate:
# admix_r <- do.call(rbind,
#                    lapply(1:5, function(i)
#                      read.table(paste0("results/NGSAdmix/", PREFIX, "/recomb_", i,"_1percent/K2.qopt")) %>%
#                        bind_cols(IDs, .) %>%
#                        left_join(., meta.ind, by = "ID") %>%
#                        arrange(., popN) %>%
#                        arrange(., zea) %>%
#                        arrange(., symp_allo) %>%
#                        filter(est_coverage > .05) %>%
#                        dplyr::mutate(., recomb_bin = paste0('recomb_', i)) %>%
#                        dplyr::mutate(., # label 'mexicana' ancestry least common in allopatric maize
#                                      mexicana_ancestry = sapply(1:nrow(.), function(j) ifelse(mean(filter(., group == "allopatric_maize")$V1) < .5, 
#                                                                                               .$V1[j], 
#                                                                                               .$V2[j])))
#                    ))
# # plot across recombination bins:
# admix_r %>%
#   dplyr::mutate(est_coverage = ifelse(est_coverage > 2, 2, est_coverage)) %>%
#   filter(est_coverage > .05) %>%
#   ggplot(aes(x = recomb_bin, y = mexicana_ancestry, 
#              color = LOCALITY, size = est_coverage)) +
#   geom_point() +
#   facet_wrap(~group) +
#   ggtitle("individual K=2 NGSAdmix mex-like ancestry estimate by r bin") +
#   xlab("low to high recombination bins (quintiles of 10kb windows)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave("plots/ind_mexicana-like-ancestry_by_recomb_bin.png",
#        height = 8, width = 8,
#        units = "in", device = "png")
# 
# # plot across recombination bins:
# admix_r %>%
#   dplyr::mutate(est_coverage = ifelse(est_coverage > 2, 2, est_coverage)) %>%
#   filter(est_coverage > 0.05) %>%
#   ggplot(aes(x = recomb_bin, y = mexicana_ancestry, 
#              color = LOCALITY, group = ID)) +
#   geom_line() +
#   #geom_point(aes(size = est_coverage)) +
#   facet_wrap(~group) +
#   ggtitle("individual K=2 NGSAdmix mex-like ancestry estimate by r bin") +
#   xlab("low to high recombination bins (quintiles of 10kb windows)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave("plots/ind_mexicana-like-ancestry_by_recomb_bin_lines.png",
#        height = 8, width = 8,
#        units = "in", device = "png")
# 
# # does the strength of selection against mexicana ancestry seem stronger 
# # for lower elevation populations?
# admix_r %>%
#   dplyr::mutate(est_coverage = ifelse(est_coverage > 2, 2, est_coverage)) %>%
#   filter(est_coverage > 0.05) %>%
#   ggplot(aes(x = recomb_bin, y = mexicana_ancestry, 
#              color = ELEVATION, group = ID)) +
#   geom_line() +
#   #geom_point(aes(size = est_coverage)) +
#   facet_wrap(~group) +
#   ggtitle("individual K=2 NGSAdmix mex-like ancestry estimate by r bin") +
#   xlab("low to high recombination bins (quintiles of 10kb windows)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave("plots/ind_mexicana-like-ancestry_by_recomb_bin_lines_colorByElevation.png",
#        height = 8, width = 8,
#        units = "in", device = "png")
# 
# # separate pattern by population -- seems to hold across most populations individually
# admix_r %>%
#   filter(., est_coverage > 0.05) %>%
#   ggplot(aes(fill = recomb_bin, x = reorder(LOCALITY, ELEVATION), y = mexicana_ancestry)) +
#   geom_boxplot() +
#   facet_wrap(~group) +
#   ggtitle("NGSadmix population mex-like ancestry estimates by r bin (boxplot of ind admixture proportions)") +
#   xlab("low to high recombination bins (quintiles of 10kb windows)") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# ggsave("plots/pop_mexicana-like-ancestry_by_recomb_bin_pops_by_locality.png",
#        height = 10, width = 12,
#        units = "in", device = "png")
# 
# admix_r %>%
#   filter(., est_coverage > 0.05) %>%
#   ggplot(aes(fill = recomb_bin, y = mexicana_ancestry)) +
#   geom_boxplot() +
#   facet_wrap(~group) +
#   ggtitle("NGSadmix mex-like ancestry estimates by r bin (boxplot of ind admixture proportions)") +
#   xlab("low to high recombination bins (quintiles of 10kb windows)") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# ggsave("plots/pop_mexicana-like-ancestry_by_recomb_bin_boxplot.png",
#        height = 6, width = 8,
#        units = "in", device = "png")
# # hmm .. looks like everyone has more mexicana ancestry in high r windows
# # including mexicana
# 
# # remake plot of ancestry ~ elevation for different recombination bins
# admix_r %>%
#   filter(., symp_allo == "sympatric") %>%
#   ggplot(., aes(x = ELEVATION, y = mexicana_ancestry, color = zea, size = est_coverage)) +
#   geom_point() +
#   ylab("mexicana ancestry") +
#   scale_colour_manual(values = colors_maize2mex[2:3]) +
#   geom_smooth(method = "lm", aes(group = zea)) +
#   ggtitle("Higher mexicana ancestry at higher elevations") +
#   facet_wrap(~recomb_bin)
# ggsave("plots/lm_predict_NGSadmix_proportion_mexicana-like_by_elevation_andByRecombBin.png", 
#        device = "png", 
#        width = 12, height = 8, units = "in",
#        dpi = 200)
# 
# # color this plot by LOCALITY
# admix_r %>%
#   filter(., symp_allo == "sympatric") %>%
#   ggplot(., aes(x = ELEVATION, y = mexicana_ancestry, color = LOCALITY, size = est_coverage)) +
#   geom_point() +
#   ylab("mexicana ancestry") +
#   geom_smooth(method = "lm", aes(group = zea), color = "black") +
#   ggtitle("Higher mexicana ancestry at higher elevations") +
#   facet_wrap(~recomb_bin)
# ggsave("plots/lm_predict_NGSadmix_proportion_mexicana-like_by_elevation_andByRecombBin_colorByPop.png", 
#        device = "png", 
#        width = 12, height = 8, units = "in",
#        dpi = 200)
# 
# 
# 
# 
# 
# # now look at Fst for the different r bins
# freqs_f <- do.call(rbind,
#                    lapply(1:5, function(i)
#                      read.table(paste0("results/NGSAdmix/", PREFIX, "/recomb_", i,"_1percent/K2.fopt.gz")) %>%
#                        dplyr::mutate(., het1 = 2*V1*(1-V1)) %>%
#                        dplyr::mutate(., het2 = 2*V2*(1-V2)) %>%
#                        dplyr::mutate(., pTot = (V1+V2)/2) %>%
#                        dplyr::mutate(., hetTot = 2*pTot*(1-pTot)) %>%
#                        dplyr::mutate(., piBetween = V1*(1-V2) + (1-V1)*V2) %>%
#                        dplyr::mutate(., recomb_bin = paste0("recomb_", i))
#                    ))
# 
# maize_is_anc1 <- admix_r %>%
#   filter(., group == "allopatric_maize") %>%
#   group_by(recomb_bin) %>%
#   dplyr::summarise(isTrue = mean(V1) > .5)
#                                                                          
# freqs_f %>%
#   dplyr::group_by(recomb_bin) %>%
#   dplyr::summarise(., Fst = 1 - (mean(het1) + mean(het2))/2/mean(hetTot),
#                    het1 = mean(het1),
#                    het2 = mean(het2),
#                    hetTot = mean(hetTot),
#                    piBetween = mean(piBetween)) %>%
#   left_join(., maize_is_anc1, by = "recomb_bin") %>%
#   dplyr::mutate(hetMaizeAnc = ifelse(maize_is_anc1$isTrue, het1, het2)) %>%
#   dplyr::mutate(hetMexicanaAnc = ifelse(maize_is_anc1$isTrue, het2, het1)) %>%
#   gather(., "stat", "value", c("hetMaizeAnc", "hetMexicanaAnc", "piBetween", "Fst", "hetTot")) %>%
#   ggplot(aes(x = recomb_bin, y = value, group = stat, color = stat)) +
#   geom_point() +
#   geom_line() +
#   facet_wrap(~stat=="Fst", scales = "free") +
#   ggtitle("pi and Fst for ancestry-estimated-allele-freqs from NGSadmix K=2 clusters") +
#   xlab("low to high recomb bins (quintiles of 10kb windows)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave("plots/fst_by_recomb_bin.png",
#        height = 6, width = 12,
#        units = "in", device = "png")
# 
# # make plots for pruned every 100th SNP:
# K = 2
# PREFIX = "pass2_alloMAIZE"
# #K = 3
# #PREFIX = "pass2_alloMAIZE_PalmarChico"
# IDs <- data.frame(ID = read.table(paste0("../samples/", PREFIX, "_IDs.list"), header = F, stringsAsFactors = F)$V1, stringsAsFactors = F)
# ancestries <- c("maize", "mexicana", "parviglumis")[1:K]
# 
# 
# file_prefix_by100 <- paste0("results/NGSAdmix/", PREFIX, "/prunedBy100/K", K)
# admix_by100 <- read.table(paste0(file_prefix_by100, ".qopt"))
# colnames(admix_by100) <- paste0("anc", 1:K) #c("anc1", "anc2")
# d_by100 <- bind_cols(IDs, admix_by100)  %>%
#   left_join(., meta, by = "ID") %>%
#   arrange(., popN) %>%
#   arrange(., zea) %>%
#   arrange(., symp_allo) %>%
#   filter(., est_coverage >= .1) %>% # only include individuals with at least .05x coverage
#   .[with(., order(group, ELEVATION)), ]
# 
# # assign mexicana ancestry
# which_anc_by100 <- data.frame(ancestry = colnames(admix_by100),
#                         ancestry_label = 
#                           sapply(colnames(admix_by100), function(x) 
#                             names(which.max(tapply(d_by100[ , x], d_by100$zea, mean)))),
#                         stringsAsFactors = F)
# 
# anc_by100 = d_by100 %>% 
#   tidyr::gather(., "ancestry", "p", colnames(admix_by100)) %>%
#   left_join(., which_anc_by100, by = "ancestry") %>%
#   dplyr::select(-ancestry) %>%
#   tidyr::spread(., ancestry_label, p)
# 
# # simple admixture plot:
# anc_by100 %>% 
#   arrange(., popN) %>%
#   tidyr::gather(., "ancestry", "p", ancestries) %>%
#   ggplot(., aes(fill=ancestry, y=p, x=ID)) +
#   geom_bar(stat = "identity", position = "fill") + 
#   facet_wrap(~group)
# 
# # plot linear models of ancestry over elevation
# lmZeaElev_by100 <- anc_by100 %>%
#   filter(., symp_allo == "sympatric") %>%
#   lm(data = ., mexicana ~ zea*ELEVATION)
# 
# anc_by100 %>%
#   filter(., symp_allo == "sympatric") %>%
#   ggplot(., aes(x = ELEVATION, y = mexicana, color = LOCALITY, size = est_coverage, shape = zea)) +
#   geom_point() +
#   ylab("Proportion mexicana ancestry") +
#   xlab("Elevation (m)") +
#   geom_abline(intercept = lmZeaElev_by100$coefficients["(Intercept)"] + lmZeaElev_by100$coefficients["zeamexicana"], 
#               slope = lmZeaElev_by100$coefficients["ELEVATION"] + lmZeaElev_by100$coefficients["zeamexicana:ELEVATION"]) +
#               #color = colors_maize2mex[3]) +
#   geom_abline(intercept = lmZeaElev_by100$coefficients["(Intercept)"], 
#               slope = lmZeaElev_by100$coefficients["ELEVATION"]) +
#               #color = colors_maize2mex[2]) +
#   ggtitle("Clines in mexicana ancestry across elevation") +
#   #theme(legend.position="bottom") +
#   labs(color = "Location", shape = "Subspecies", size = "Coverage") +
#   guides(color = FALSE) +
#   theme_classic()
# ggsave(paste0("plots/lm_predict_NGSadmix_proportion_mexicana-like_by_elevation_colored_by_pop_prunedBy100_", PREFIX, "_K", K, ".pdf"), 
#        device = "pdf", 
#        width = 5, height = 4, units = "in")
# ggsave(paste0("../../hilo_manuscript/figures/lm_predict_NGSadmix_proportion_mexicana-like_by_elevation_colored_by_pop_prunedBy100_", PREFIX, "_K", K, ".pdf"), 
#        device = "pdf", 
#        width = 5, height = 4, units = "in")
# 

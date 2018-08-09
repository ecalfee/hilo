library(scales)
library(dplyr)
library(tidyr)
library(ggplot2)

# ID and population labels (created by addMetadata2HiloIDs.R):
pass1 <- read.table("../data/pass1_ids.txt", stringsAsFactors = F, header = T, sep = "\t")

# first visualize basic depth calculation in a histogram:
png('../plots/pass1_hist_coverage_passQC.png', width = 800, height = 600)
hist(pass1$est_coverage, breaks = 20, 
     main = "pass1 est. coverage post filtering", 
     xlab = "(# reads pass QC * 150bp) / 2.1Gb")
abline(v = median(pass1$est_coverage), col = "blue")
abline(v = mean(pass1$est_coverage), col = "red")
dev.off()

# number of pass1 individuals with higher coverage
sum(pass1$est_coverage>.1)
sum(pass1$est_coverage>.25)
sum(pass1$est_coverage>.5)


# plot mapping statistics across hilo: % mapped; % mapped > mapQ30 by mex/maize:
# add in alignment metrics: mapping quality and sequencing coverage
raw_metrics_no_header <- read.table("../data/filtered_bam/pass1.all.metrics.raw.Nreads", 
                                    stringsAsFactors = F, header = F, skip = 1)
con = file("../data/filtered_bam/pass1.all.metrics.raw.Nreads", "r")
metrics_head <- strsplit(readLines(con, 1), split = "\t")[[1]]
close(con)
raw_metrics0 <- raw_metrics_no_header %>% 
  unite(., "LIB", V2, V3, sep= "_") # join 2 library_unknown columns
colnames(raw_metrics0) <- metrics_head # add header
raw_metrics <- raw_metrics0 %>%
  mutate(n_reads = 2*READ_PAIRS_EXAMINED + UNPAIRED_READS_EXAMINED) %>%
  mutate(prop_unmapped = UNMAPPED_READS/n_reads) %>%
  mutate(n_dedup = 2*(READ_PAIRS_EXAMINED - READ_PAIR_DUPLICATES) +
           UNPAIRED_READS_EXAMINED - UNPAIRED_READ_DUPLICATES) %>%
  mutate(n = as.integer(substr(StudyID, 5, 10))) %>%
  left_join(., pass1, by = "n") %>%
  mutate(n_filtOutQ30 = n_dedup-num_read_pass) %>%
  mutate(prop_filtOutQ30 = n_filtOutQ30/n_dedup) %>%
  mutate(prop_filtTot = (n_reads - num_read_pass)/n_reads) %>%
  mutate(n_dropped = n_reads-n_dedup+n_filtOutQ30)

# average total depth by group
raw_metrics %>%
  group_by(paste(zea, symp_allo, sep = "_")) %>%
  summarise(., total = sum(est_coverage))
raw_metrics %>%
  group_by(paste(zea, symp_allo, sep = "_")) %>%
  summarise(., total = sum(depthRegions))
# total coverage from allopatric mexicana is <6x
# only 12 of the 28 individuals have individual coverage >.2
raw_metrics %>%
  group_by(paste(zea, symp_allo, sep = "_")) %>%
  count(., est_coverage > .2)

# Below I'm analysing coverage data pulled from ANGSD
# from output files of calcDepthCovSites.sh ()
# S for sites; these files are split by chromosome
# chr1 excluded for now (Did not finish running)
depthS_nochr1 <- Reduce('+', 
                        lapply(2:10, function(chr)
                          read.csv(paste0("../data/depthCov/pass1/all.positions.chr",
                                          chr, ".depthSample"),
                                   header = F, sep = "\t")))
depthS20_nochr1 <- Reduce('+', 
                          lapply(2:10, function(chr)
                            read.csv(paste0("../data/depthCov/pass1/all.positions.chr",
                                            chr, ".Q20.depthSample"),
                                     header = F, sep = "\t")))

# and calcDepthCovRegions.sh (regions used to calc Fst)
# F for Fst
depthF <- read.csv("../data/depthCov/pass1/N1000.L100.regions.depthSample",
                   header = F, sep = "\t")
depthF20 <- read.csv("../data/depthCov/pass1/N1000.L100.regions.Q20.depthSample",
                     header = F, sep = "\t")

# function to calculate mean depth from raw depth dataframe
meanDepth = function(dep){
  dep$V102 <- NULL
  colnames(dep)<- paste0("n", 0:100)
  meanDepth <- apply(dep[ , paste0("n", 0:100)], 1, 
                         function(r) r %*% 0:100/sum(r))
  return(meanDepth)
}
# calculate and write to file mean depth across the Fst regions from depthF
mean_depthF = meanDepth(depthF)
plot(pass1$est_coverage, mean_depthF, main = "Depth est. regions vs. # reads") # perfectly correlated
write.table(mean_depthF, "../data/depthCov/pass1/N1000.L100.regions_meanDepth.txt", 
            col.names = T, row.names = T, quote = F)


# write a new metrics file with depthF and all other filtering measures too
raw_metrics$depthRegions = mean_depthF
write.table(raw_metrics, "../data/filtered_bam/pass1.all.metrics.calcs", 
            col.names = T, row.names = F, quote = F)




# plots for coverage and filtering effects
png(paste("../plots/pass1_total_reads_dropped_by_group.png"), # saves plot as pin in ../plots/
    height = 6, width = 6, units = "in", res = 150)
qplot(prop_filtTot, data=raw_metrics, geom="density", fill=paste(zea, symp_allo, sep = "_"), alpha=I(.5), 
      main="Prop. reads filtered total by group", xlab="Prop. reads dropped", 
      ylab="Density") +
  theme(legend.title=element_blank())
dev.off()
# same but in scatterplot across sequencing coverage
ggplot(data = raw_metrics) + geom_point(aes(x = n_reads, 
                                            y = prop_filtTot,
                                            color = paste(zea, symp_allo, sep = "_"))) +
  theme(legend.title=element_blank())
png(paste("../plots/pass1_reads_unmapped_by_group.png"), # saves plot as pin in ../plots/
    height = 6, width = 6, units = "in", res = 150)
qplot(prop_unmapped, data=raw_metrics, geom="density", fill=paste(zea, symp_allo, sep = "_"), alpha=I(.5), 
      main="Prop. unmapped reads by group", xlab="Prop. reads unmapped", 
      ylab="Density") +
  theme(legend.title=element_blank())
dev.off()
ggplot(data = raw_metrics) + geom_point(aes(x = n_reads, 
                                            y = prop_unmapped,
                                            color = paste(zea, symp_allo, sep = "_"))) +
  theme(legend.title=element_blank())
# of reads that mapped and are not duplicates -- what proportion are mapQ <30?
png(paste("../plots/pass1_prop_reads_mapped_but_belowQ30_by_group.png"), # saves plot as pin in ../plots/
    height = 6, width = 6, units = "in", res = 150)
qplot(prop_filtOutQ30, data=raw_metrics, geom="density", fill=paste(zea, symp_allo, sep = "_"), alpha=I(.5), 
      main="Prop. reads filtered by mapQ<30 by group", xlab="Prop. reads < Q30", 
      ylab="Density") +
  theme(legend.title=element_blank())
dev.off()
# strangely, teosinte allopatric has higher percent duplications
png(paste("../plots/pass1_prop_duplication_by_group.png"), # saves plot as pin in ../plots/
    height = 6, width = 6, units = "in", res = 150)
qplot(PERCENT_DUPLICATION, data=raw_metrics, geom="density", fill=paste(zea, symp_allo, sep = "_"), alpha=I(.5), 
      main="Prop. reads duplications by group", xlab="Prop. reads duplicates", 
      ylab="Density") +
  theme(legend.title=element_blank())
dev.off()

# unclear what the effect of BAQ (adjusted base quality) filtering is .. could try to count for a set of SNPs
# Jeff thinks a larger issue is getting total coverage across individuals, groups and entire sample for SNPs included (e.g. in PCA)

# function to plot mean depth
plot_depth_mean = function(depData, ids, label){ # takes in depth dataframe, matching individual ID's (e.g. pass1), and a label for plots
  dep = ids 
  dep$meanDepth = meanDepth(depData)
  plot2 = qplot(meanDepth, data=dep, geom="density", 
        fill=paste(zea, symp_allo, sep = "_"), alpha=I(.5)) +
    labs(x="mean depth", 
         y="Density",
         title=label,
         subtitle="Depth mapQ>30 by group") +
    theme(legend.title=element_blank())
  print(plot2) # ggplots must be called explicitly to print
ggsave(paste0("../plots/hist of individuals -- mean", label, ".png"), plot = plot2, device = png())
  }
# plot for all
depths_list = list(depthF, depthF20, depthS_nochr1, depthS20_nochr1)
depths_labels = c("depth Fst regions", "depth Fst regions Q>20",
                  "depth var sites", "depth var sites Q>20")
for (i in 1:4) plot_depth_mean(dep = depths_list[[i]], ids = pass1, label = depths_labels[i])

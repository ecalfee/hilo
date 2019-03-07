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
sum(pass1$est_coverage>.05)
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
#depthF <- read.csv("../data/depthCov/pass1/N1000.L100.regions.depthSample",
#                   header = F, sep = "\t")
#depthF20 <- read.csv("../data/depthCov/pass1/N1000.L100.regions.Q20.depthSample",
#                     header = F, sep = "\t")
depthF <- read.csv("../variant_sites/results/depthCov/N1000.L100.regions/hilo_alloMAIZE_MAIZE4LOW.depthSample",
                   header = F, sep = "\t")
depthF20 <- read.csv("../variant_sites/results/depthCov/N1000.L100.regions/hilo_alloMAIZE_MAIZE4LOW.Q20.depthSample",
                     header = F, sep = "\t")
# global depth
depthG <- read.csv("../variant_sites/results/depthCov/N1000.L100.regions/hilo_alloMAIZE_MAIZE4LOW.depthGlobal",
                   header = F, sep = "\t")
depthG20 <- read.csv("../variant_sites/results/depthCov/N1000.L100.regions/hilo_alloMAIZE_MAIZE4LOW.Q20.depthGlobal",
                     header = F, sep = "\t")


# function to calculate mean depth from raw depth dataframe
meanDepth = function(dep, max = 100){
  dep[, ncol(dep)] <- NULL
  colnames(dep)<- paste0("n", 0:max)
  meanDepth <- apply(dep[ , paste0("n", 0:max)], 1, 
                         function(r) r %*% (0:max)/sum(r))
  return(meanDepth)
}
# what is the mean global depth?
meanDepth(depthG20, max = 10000)
meanDepth(depthG, max = 10000)
sum(meanDepth(depthF))
sum(meanDepth(depthF20))


# calculate and write to file mean depth across the Fst regions from depthF
mean_depthF = meanDepth(depthF)
plot(pass1$est_coverage, mean_depthF, main = "Depth est. regions vs. # reads") # perfectly correlated

# average total depth by group
raw_metrics %>%
  group_by(paste(zea, symp_allo, sep = "_")) %>%
  summarise(., total = sum(est_coverage))
raw_metrics %>%
  mutate(mean_depthF = mean_depthF) %>%
  group_by(paste(zea, symp_allo, sep = "_")) %>%
  summarise(., total = sum(mean_depthF))
# total coverage from allopatric mexicana is <6x
# only 12 of the 28 individuals have individual coverage >.2
raw_metrics %>%
  group_by(paste(zea, symp_allo, sep = "_")) %>%
  count(., est_coverage > .2)
raw_metrics %>%
  mutate(mean_depthF = mean_depthF) %>%
  group_by(paste(zea, symp_allo, sep = "_")) %>%
  count(., mean_depthF > .2)


write.table(mean_depthF, "../data/depthCov/pass1/N1000.L100.regions_meanDepth.txt", 
            col.names = T, row.names = T, sep = "\t", quote = F)


# write a new metrics file with depthF and all other filtering measures too
raw_metrics$depthRegions = mean_depthF
write.table(raw_metrics, "../data/filtered_bam/pass1.all.metrics.calcs", 
            col.names = T, row.names = F, sep = "\t", quote = F)




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

# get total depth of coverage distribution for hilo ind's
# combined wiht maize 4 low (16 higher cov. ind's) from Anne
depthT <- t(read.csv("../data/depthCov/merged_pass1_all_alloMaize4Low_16/N1000.L100.regions.Q20.depthGlobal",
                   header = F, sep = "\t"))[1:5001]
depthBin <- 0:5000
depthT[1:5] # no spots with zero coverage but 867 with 1x coverage?
depthTvals <- unlist(mapply(rep, times = depthT, x = depthBin))
hist(depthTvals)
summary(depthTvals)
sd(depthTvals)
mean(depthTvals)
maxDepth = mean(depthTvals) + 3*sd(depthTvals)
# 1020
# how many sites does this cutoff? ~ 1.2%
sum(depthTvals>maxDepth)/length(depthTvals)
# depth per sample (just hilo, not allopatric maize)
depthS <- read.csv("../data/depthCov/merged_pass1_all_alloMaize4Low_16/N1000.L100.regions.Q20.depthSample",
                     header = F, sep = "\t")[, 1:5001]
depthShilo <- depthS[1:196,]
apply(depthShilo, 2, sum)[1:50] # hilo totals
depthSvalshilo <- lapply(1:nrow(depthShilo), function(r) unlist(mapply(rep, times = depthShilo[r,], x = depthBin)))
depthSvalsAllhilo <- unlist(depthSvalshilo)
hist(depthSvalsAllhilo)
summary(depthSvalsAllhilo)
mean(depthSvalsAllhilo) + 3*sd(depthSvalsAllhilo)
summary(sapply(depthSvalshilo, mean) + 3*sapply(depthSvalshilo, sd))
# 8 is an upper cutoff for 3*sd for an individual
summary(sapply(depthSvalshilo, function(x) sum(x>8)/length(x)))
# for the highest coverage individual this is ~1.4% of sites
summary(sapply(depthSvalshilo, mean) + 6*sapply(depthSvalshilo, sd))
mean(depthSvalsAllhilo) + 10*sd(depthSvalsAllhilo)
# 12 is a very conservative cutoff
mean(depthSvalsAllhilo) + 7*sd(depthSvalsAllhilo)

# only allopatric maize
depthSmaize <- depthS[197:212,]
apply(depthSmaize, 2, sum)[1:50] # maize totals
depthSvalsmaize <- lapply(1:nrow(depthSmaize), function(r) unlist(mapply(rep, times = depthSmaize[r,], x = depthBin)))
depthSvalsAllmaize <- unlist(depthSvalsmaize)
hist(depthSvalsAllmaize)
summary(depthSvalsAllmaize)
mean(depthSvalsAllmaize) + 3*sd(depthSvalsAllmaize)
summary(sapply(depthSvalsmaize, mean) + 3*sapply(depthSvalsmaize, sd))
# 63 is an upper cutoff for 3*sd for an individual
summary(sapply(depthSvalsmaize, function(x) sum(x>63)/length(x)))
# for the highest coverage individual, this is ~1% of sites
mean(depthSvalsAllmaize) + 4*sd(depthSvalsAllmaize)
mean(depthSvalsAllmaize) + 7*sd(depthSvalsAllmaize)
mean(depthSvalsAllmaize) + 10*sd(depthSvalsAllmaize)
# 158 is a very conservative cutoff allowing most sites to pass filtering

# what percent of sites do we expect to have zero coverage for an individual?
pdataHilo = 1 - sum(depthSvalsAllhilo == 0)/length(depthSvalsAllhilo)
pdataMaize = 1 - sum(depthSvalsAllmaize == 0)/length(depthSvalsAllmaize)
# what's the probability of at least 6 ind's (>4) with data 
# out of 34 hilo ind's in mex reference panel?
pbinom(q = 5, size = 34, prob = pdataHilo, lower.tail = F)
# 12 individuals in maize is fairly similar probability
pbinom(q = 11, size = 16, prob = pdataMaize, lower.tail = F)
# less strict cutoffs that about 86% sites should pass
# in either pop are minInd_mex = 5 and minInd_maize=11
pbinom(q = 4, size = 34, prob = pdataHilo, lower.tail = F)
pbinom(q = 10, size = 16, prob = pdataMaize, lower.tail = F)

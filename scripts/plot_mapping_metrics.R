library(scales)
library(dplyr)
library(tidyr)
library(ggplot2)

# ID and population labels:
hilo <- read.table("../data/HILO_IDs_cov_pass1.csv", stringsAsFactors = F, header = T)
pass1 <- hilo[!is.na(hilo$num_read_pass),] # only include individuals with some pass1 data


# plot mapping statistics across hilo: % mapped; % mapped > mapQ30 by mex/maize:
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
  left_join(., hilo, by = "n") %>%
  mutate(n_filtOutQ30 = n_dedup-num_read_pass) %>%
  mutate(prop_filtOutQ30 = n_filtOutQ30/n_dedup) %>%
  mutate(prop_filtTot = (n_reads - num_read_pass)/n_reads) %>%
  mutate(n_dropped = n_reads-n_dedup+n_filtOutQ30)


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

# function to tidy the data
tidy_depth = function(dep){
  dep$V102 <- NULL
  colnames(dep)<- paste0("n", 0:100)
  dep <- cbind(pass1, dep) # add hilo id's & basic metadata
  dep$meanDepth <- apply(dep[ , paste0("n", 0:100)], 1, 
                         function(r) r %*% 0:100/sum(r))
  dep2 = gather(dep, "depth", "n_sites", paste0("n", 0:100))
  dep3 = separate(dep2, depth, c("n", "site_depth"), sep = 1) %>%
    select(., -n) %>%
    mutate(., site_depth = as.integer(site_depth))
}
group_by(dep3, c(zea, symp_allo)) %>% summarise(., total_by_grp = sum(n_sites))
ggplot(data = dep3, aes(site_depth, n_sites)) + geom_point(aes(color = zea))

# function to plot mean depth
plot_depth_mean = function(dep, label){ # takes in depth dataframe and a label for plots
  dep$V102 <- NULL
  colnames(dep)<- paste0("n", 0:100)
  dep <- cbind(pass1, dep) # add hilo id's & basic metadata
  dep$meanDepth <- apply(dep[ , paste0("n", 0:100)], 1, 
                         function(r) r %*% 0:100/sum(r))
  
  plot(x = 0:100, y = seq(0, 100000, length.out = 101), 
       col = NULL, xlab = "depth at site",
       ylab = "number of sites",
       main = label,
       sub = "depth of reads mapped >30",
       xlim = c(0,20))
  for (i in 1:nrow(dep)){
    points(0:100, dep[i, paste0("n", 0:100)], type = "l", 
           col = ifelse(dep[i, "zea"]=="maize", "yellow", "blue"))
  }

  
  plot2 = qplot(meanDepth, data=dep, geom="density", 
        fill=paste(zea, symp_allo, sep = "_"), alpha=I(.5)) +
    labs(x="mean depth", 
         y="Density",
         title=label,
         subtitle="Depth mapQ<30 by group") +
    theme(legend.title=element_blank())
  print(plot2) # ggplots must be called explicitly to print
}
# plot for all
depths_list = list(depthF, depthF20, depthS_nochr1, depthS20_nochr1)
depths_labels = c("depth Fst regions", "depth Fst regions Q>20",
                  "depth var sites", "depth var sites Q>20")
for (i in 1:4) plot_depth_mean(dep = depths_list[[i]], label = depths_labels[i])


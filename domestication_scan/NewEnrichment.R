library(GenomicRanges)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(seqinr)

# script to calculate bp overlap between raisd putative domestication regions
# and 'introgression deserts' in sympatric maize and mexicana 
# working directory is hilo/

# input files
raisd_input <- "domestication_scan/results/RAiSD_Output.Raisd"
rmap_file <- "data/linkage_map/ogut_fifthcM_map_agpv4_EXTENDED.txt"
maize_bed <- "local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/maize.combined.anc.bed"
mexicana_bed <- "local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne10000_yesBoot/anc/mexicana.combined.anc.bed"
# output files
ref_maize <- "data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
fa_out <- "domestication_scan/results/out.fa"
raisd_bed <- "domestication_scan/results/raisdHits.bed" 

  
#load in Raisd results
raisd<-fread(raisd_input,data.table = FALSE)
colnames(raisd)<-c('seqnames','testsnp','start','end','a','b','c','score')

raisd$seqnames<-paste('chr',raisd$seqnames,sep="")
raisdraw<-raisd #for bootstrap test

#pull in the genetic map
map<-fread(rmap_file,data.table=FALSE)
colnames(map)<-c('SNP','marker','pos_cM', 'chr', 'pos_bp')

#impute recomb rate for raisd snps using linear interpolation:
raisdNew <-do.call(rbind,
                           lapply(1:10,function(cr){
  raisd_subset<-subset(raisd,seqnames==paste('chr',cr,sep=""))
  map_subset<-subset(map,chr==cr)
  raisd_subset$pos_cm <- approx(x = map_subset$pos_bp,
                                y = map_subset$pos_cM,
                                xout = raisd_subset$testsnp,
                                method = "linear")$y
   return(raisd_subset)
}))

#raisdNew0<-bind_rows(subsets)
#raisdNew0<-raisdNew0[-which(abs(raisdNew0$pos_cm)>500),] #remove recombination outliers for pairing (around 1%)

#now chunk recombs into even deciles and sample windows evenly:
bin<-as.numeric(cut_number(raisdNew$pos_cm,10))
raisdNew$bin<-bin

# plot: these bins are cM position deciles across genome (not cM/bp recombination rates)
test_bins <- data.frame(bin = bin[c(T, rep(F, 999))], 
                        bp_pos = raisdNew$testsnp[c(T, rep(F, 999))],
                        chr = raisdNew$seqnames[c(T, rep(F, 999))])
ggplot(test_bins, aes(x = bp_pos, y = bin, color = chr)) +
  geom_point() +
  theme_classic() +
  facet_grid(chr~.)

#calculate window width for matching
raisdNew$width<-raisdNew$end-raisdNew$start

#Now past the heavy lifting on raisd's end (make a save and load for dev)


#do the actual subsetting here by score (raisd score)
#This sort is memory intensive, doing it in two steps
#First remove the 0 scores first (we'll never select anyhow)
raisdNew2<-raisdNew[-which(raisdNew$score==0),]

#Then do the sort
raisdSort<-raisdNew2[order(raisdNew2$score,decreasing=TRUE),]
#Rough attempt to get top X%
raisdHits<-raisdSort[1:round(nrow(raisdNew)*0.001),] #testing because the final length will depend on how much overlap is here

grRaisd <- makeGRangesFromDataFrame(raisdHits, keep.extra.columns=TRUE) %>%
  reduce()

#write a bed file of grRaisd and remove the ones that are mostly Ns
rdf<-data.frame(grRaisd)
rdf$seqnames<-sapply(rdf$seqnames,substr,4,6) # chr2 -> 2

write.table(rdf[,c(1:3)],file=raisd_bed,
            quote = FALSE,row.names = FALSE,col.names = FALSE,sep = '\t')

#now run bedtools:
system(paste("bedtools getfasta -fi", ref_maize, "-fo", fa_out, "-bed", raisd_bed))
       
#now check to see which of the outliers are driven by N's

#read in the bed file:
hits<-fread(raisd_bed,data.table = FALSE)

#read in the sequences of the regions:
fas<-read.fasta(fa_out)
totals<-sapply(fas,function(x) length(which(x == 'n')))
percents<-sapply(fas,function(x) length(which(x == 'n'))/length(x))
#Use for weighting bootstraps
kps<-rdf[which(percents<0.5),] # keep only raisd hits that have less than 50% n's (missing data in reference genome)
kps$seqnames<-sapply(kps$seqnames,function(x) paste('chr',x,sep=""))
perc<-0.7
#remake the genomicranges file:
grRaisd_kps <- makeGRangesFromDataFrame(kps, keep.extra.columns=TRUE)


#load in introgression scan data:
# sympatric maize data (mean mexicana ancestry at each position)
maizraw<-fread(maize_bed,data.table = FALSE)
colnames(maizraw)<-c('seqnames','start','end','testsnp','score')
maizraw$seqnames<-paste('chr',maizraw$seqnames,sep="")

#find the break for candidates
twop<-round(0.02*dim(maizraw)[1]) # changed from .005 to .02 for 2% (not 0.5%)

# lowest 2% mexicana introgression in maize
cutoff<-maizraw$score[order(maizraw$score)[twop]]
# table(maiz$score<cutoff)/nrow(maiz)*100 # 2%
#subset by score (mean mex ancestry)
maiz<-maizraw[which(maizraw$score<cutoff),]

# sympatric mexicana data (mean mexicana ancestry at each position)
mexraw<-fread(mexicana_bed,data.table = FALSE)
colnames(mexraw)<-c('seqnames','start','end','testsnp','score')
mexraw$seqnames<-paste('chr',mexraw$seqnames,sep="")
#find the X% break
twop<-round(0.02*dim(mexraw)[1]) # changed from 0.06 to 0.02 for 2%

# highest 2% mexicana ancestry (low introgression) in mexicana
cutoff<-mexraw$score[order(mexraw$score,decreasing = TRUE)[twop]]
#subset by score (mean mex ancestry)
# table(mexraw$score>cutoff)/nrow(mexraw)*100 # 2%
mex<-mexraw[which(mexraw$score>cutoff),]


#cM spline for maiz/mex

MaizNew <- do.call(rbind,
                   lapply(1:10,function(cr){
  maiz_subset<-subset(maiz,seqnames==paste('chr',cr,sep=""))
  map_subset<-subset(map,chr==cr)
  maiz_subset$pos_cm <- approx(x = map_subset$pos_bp,
                               y = map_subset$pos_cM,
                               xout = maiz_subset$testsnp,
                               method = "linear")$y
  return(maiz_subset)
}))

MaizNew$width<-MaizNew$end-MaizNew$start

grMaiz <- makeGRangesFromDataFrame(MaizNew, keep.extra.columns=TRUE)


#be aware that based on filtering you may not have hits on a chromosome and that will trip up things downstream
MexNew <- do.call(rbind, 
                  lapply(c(1:10),function(cr){
  mex_subset<-subset(mex,seqnames==paste('chr',cr,sep=""))
  map_subset<-subset(map,chr==cr)
  mex_subset$pos_cm <- approx(x = map_subset$pos_bp,
                              y = map_subset$pos_cM,
                              xout = mex_subset$testsnp,
                              method = "linear")$y
  return(mex_subset)
}))

MexNew$width<-MexNew$end-MexNew$start

grMex <- makeGRangesFromDataFrame(MexNew, keep.extra.columns=TRUE)

#now compare overlaps

#Get an estimate of the sizes of tests that I'm looking at to get an expectation based on percentage


hitsMex <- countOverlaps(grRaisd_kps, grMex, ignore.strand=TRUE)
#I don't think I want sum, I think I want is 0 or not
length(which(hitsMex>0))
raisdlen<-length(grRaisd_kps)
raisdlenbases<-sum(data.frame(grRaisd_kps)$width)
Mexlen<-length(grMex)
sum(data.frame( GenomicRanges::intersect(grRaisd_kps, grMex,ignore.strand=TRUE))$width)

hitsMaiz <- countOverlaps(grRaisd_kps, grMaiz, ignore.strand=TRUE)
length(which(hitsMaiz>0))
maizlen<-length(grMaiz)
sum(data.frame( GenomicRanges::intersect(grRaisd_kps, grMaiz,ignore.strand=TRUE))$width)

hitsMexMaiz <- countOverlaps(grMex, grMaiz, ignore.strand=TRUE)
length(which(hitsMexMaiz>0))
sum(data.frame( GenomicRanges::intersect(grMex, grMaiz,ignore.strand=TRUE))$width)

#get approx proportions of bins from final:
proportions<-table(raisd$bin)/nrow(raisd)


MaizNew$bin<-as.numeric(cut_number(MaizNew$pos_cm,10))
maiz<-MaizNew
maizproportions<-table(maiz$bin)

MexNew$bin<-as.numeric(cut_number(MexNew$pos_cm,10))
mex<-MexNew
mexproportions<-table(mex$bin)

#Right now I'm not doing anything with rasid/mex or mex/maize overlaps but I'm calculating them just in case
#Not certain what the expectation or biological explanation would be if there were enrichment

#now iterate the function to sample randomly and build a bootstrap distribution

#also update raisdraw seqnames once to save time:
#raisdraw$seqnames<-paste('chr',raisdraw$seqnames,sep="")

#get the window sizes of raisd:
rdf_kps<-data.frame(grRaisd_kps)


#Now the hybridization islands are held constant. I just permute the raisd window sizes onto different locations

boots<-sapply(1:1000,function(x){
  print(x)
  #load in Raisd
  raisd<-raisdraw
  colnames(raisd)<-c('seqnames','testsnp','start','end','a','b','c','score')
  
  #subset by score (raisd score)
  binvals<-sample(1:10,raisdlen,prob=proportions,replace=TRUE)
  
  lists<-lapply(1:10,function(x) sample(which(raisdNew$bin==x),10000))
  #now randomly grab a line from the respective lists:
  grabs<-sapply(binvals,function(x) sample(lists[[x]],1))

  #takes existing sizes of windows and drops them into different recomb match regions
  raisdgr<-raisdNew[grabs,]
  wid<-(rdf$width/2)*(1-perc) #correct for width and Ns in empirical
  raisdgr$start<-raisdgr$testsnp-wid
  raisdgr$end<-raisdgr$testsnp+wid
  raisd<-raisdgr
  
  grRaisd <- makeGRangesFromDataFrame(raisd, keep.extra.columns=TRUE)
  grRaisd<-reduce(grRaisd)
  
  #maiz<-maizraw
  #colnames(maiz)<-c('seqnames','start','end','testsnp','score')
  #maiz$seqnames<-paste('chr',maiz$seqnames,sep="")
  #maiz<-maiz[sample(1:nrow(maiz),nrow(MaizNew)),]
  #not included, also reshuffling Maiz/Mex (old)
  #subsets<-lapply(1:10,function(cr){
  #  maiz_subset<-subset(maiz,seqnames==paste('chr',cr,sep=""))
  #  map_subset<-subset(map,chr==cr)
  #  maiz_subset$pos_cm <-spline(x=map_subset$pos_bp,y=map_subset$pos_cM,xout=maiz_subset$testsnp, method = "hyman")$y
  #  return(maiz_subset)
  #})
  
  #MaizNew<-bind_rows(subsets)
  
  #grMaiz <- makeGRangesFromDataFrame(MaizNew, keep.extra.columns=TRUE)
  #grMaiz <- reduce(grMaiz) #this contiguous windows
  
  
  #now compare overlaps
  hitsMex <- sum(data.frame( GenomicRanges::intersect(grRaisd, grMex,ignore.strand=TRUE))$width)
  hitsMaiz <- sum(data.frame( GenomicRanges::intersect(grRaisd, grMaiz,ignore.strand=TRUE))$width)
  hitsMexMaiz <- sum(data.frame( GenomicRanges::intersect(grMex, grMaiz,ignore.strand=TRUE))$width)
  print(hitsMaiz)
  return(c(hitsMex,hitsMaiz,hitsMexMaiz))

})

apply(boots,MARGIN=1,range)

#plot the distribution of maize vs raisd hits
bootdf<-data.frame(t(boots))
colnames(bootdf)<-c('Mex','Maiz','MexMaiz')


#Now save the requisite data:
maizOverlap<-sum(data.frame( GenomicRanges::intersect(grRaisd, grMaiz,ignore.strand=TRUE))$width)
mexOverlap<-sum(data.frame( GenomicRanges::intersect(grRaisd, grMex,ignore.strand=TRUE))$width)
save(bootdf,maizOverlap,mexOverlap,mexraw,raisdraw,maizraw,file='DfForEnrichment.Rimage')



#####################################################################################################
################ Start here for plotting ############################################################
#####################################################################################################

load('DfForEnrichment.Rimage')

#sum(tapply(mexraw$start,mexraw$chr,max))
#sum(tapply(raisdraw$start,raisdraw$seqnames,max))
#both datasets are roughly 2G (genome size)

#sum(data.frame(grRaisd)$width)/sum(tapply(raisdraw$start,raisdraw$seqnames,max))
#sum(data.frame(grMaiz)$width)/sum(tapply(maizraw$start,maizraw$chr,max))
#sum(data.frame(grMex)$width)/sum(tapply(mexraw$start,mexraw$chr,max))

percentMaiz<-maizOverlap/sum(tapply(raisdraw$start,raisdraw$seqnames,max))
percentMex<-mexOverlap/sum(tapply(raisdraw$start,raisdraw$seqnames,max))

#I don't think these percents are what I really want (I think matching by total number of bases is the best)
bootdf$MexPerc<-bootdf$Mex/sum(tapply(raisdraw$start,raisdraw$seqnames,max))
bootdf$MaizPerc<-bootdf$Maiz/sum(tapply(raisdraw$start,raisdraw$seqnames,max))

png('EnrichmentPlot.png',res = 250,height=960,width=960)
ggplot(data=bootdf,aes(x=MaizPerc))+
  geom_histogram(color='blue')+
  geom_vline(xintercept = percentMaiz,linetype='dashed',color='red')+
  xlab('Number of Raisd and Hybrid Desert Overlaps')+
  ylab('Bootstrap Counts')
dev.off()


png('EnrichmentPlotMex.png',res = 250,height=960,width=960)
ggplot(data=bootdf,aes(x=MexPerc))+
  geom_histogram(color='blue')+
  geom_vline(xintercept = percentMex,linetype='dashed',color='red')+
  xlab('Number of Raisd and Hybrid Desert Overlaps')+
  ylab('Bootstrap Counts')
dev.off()



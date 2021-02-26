library(GenomicRanges)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(seqinr)
setwd('/Users/danielgates/Desktop/ErinHybridization/')

#load in Raisd
raisd<-fread('/Users/danielgates/Desktop/ErinHybridization/RAiSD_Output.Raisd',data.table = FALSE)
colnames(raisd)<-c('seqnames','testsnp','start','end','a','b','c','score')

raisd$seqnames<-paste('chr',raisd$seqnames,sep="")
raisdraw<-raisd #for bootstrap test

#pull in the genetic map
#map<-fread('ogut_fifthcM_map_agpv4_EXTENDED.txt',data.table=FALSE)
map<-fread('ogut_fifthcM_map_agpv4_INCLUDE.txt',data.table=FALSE)
colnames(map)<-c('SNP','marker','pos_cM', 'chr', 'pos_bp')

#impute recomb rate for raisd snps:
subsets<-lapply(1:10,function(cr){
  raisd_subset<-subset(raisd,seqnames==paste('chr',cr,sep=""))
  map_subset<-subset(map,chr==cr)
  raisd_subset$pos_cm <-spline(x=map_subset$pos_bp,y=map_subset$pos_cM,xout=raisd_subset$testsnp, method = "hyman")$y
  return(raisd_subset)
})

raisdNew<-bind_rows(subsets)
raisdNew<-raisdNew[-which(abs(raisdNew$pos_cm)>500),] #remove recombination outliers for pairing (around 1%)

#now chunk recombs into even deciles and sample windows evenly:
bin<-as.numeric(cut_number(raisdNew$pos_cm,10))
raisdNew$bin<-bin
raisd<-raisdNew

#calculate window width for matching
raisd$width<-raisd$end-raisd$start

#Now past the heavy lifting on raisd's end (make a save and load for dev)


#do the actual subsetting here by score (raisd score)
#This sort is memory intensive, doing it in two steps
#First remove the 0 scores first (we'll never select anyhow)
raisd2<-raisd[-which(raisd$score==0),]

#Then do the sort
raisdSort<-raisd2[order(raisd2$score,decreasing=TRUE),]
#Rough attempt to get top X%
raisd<-raisdSort[1:round(nrow(raisd)*0.001),] #testing because the final length will depend on how much overlap is here

grRaisd <- makeGRangesFromDataFrame(raisd, keep.extra.columns=TRUE)
grRaisd2<-grRaisd
grRaisd<-reduce(grRaisd)

#write a bed file of grRaisd and remove the ones that are mostly Ns
rdf<-data.frame(grRaisd)
rdf$seqnames<-sapply(rdf$seqnames,substr,4,6)

write.table(rdf[,c(1:3)],file="raisdHits.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = '\t')

#now run bedtools:
system("bedtools getfasta -fi Zea_mays.AGPv4.dna.chr.fa -fo out.fa -bed raisdHits.bed ")

#now check to see which of the outliers are driven by N's

#read in the bed file:
hits<-fread('raisdHits.bed',data.table = FALSE)

#read in the sequences of the regions:
fas<-read.fasta('out.fa')
totals<-sapply(fas,function(x) length(which(x == 'n')))
percents<-sapply(fas,function(x) length(which(x == 'n'))/length(x))
#Use for weighting bootstraps
kps<-rdf[which(percents<0.5),]
kps$seqnames<-sapply(kps$seqnames,function(x) paste('chr',x,sep=""))
perc<-0.7
#remake the genomicranges file:
grRaisd <- makeGRangesFromDataFrame(kps, keep.extra.columns=TRUE)


#load in introgression scan data:
maizraw<-fread('maize.combined.anc.bed',data.table = FALSE)
maiz<-maizraw
colnames(maiz)<-c('seqnames','start','end','testsnp','score')
maiz$seqnames<-paste('chr',maiz$seqnames,sep="")

#find the break for candidates
twop<-round(0.005*dim(maiz)[1])

#increasing:
cutoff<-maiz$score[order(maiz$score)[twop]]
#subset by score (mean mex ancestry)
maiz<-maiz[which(maiz$score<cutoff),]

#decreasing:
#cutoff<-maiz$score[order(maiz$score,decreasing = TRUE)[twop]]
#maiz<-maiz[which(maiz$score>cutoff),]

mexraw<-fread('mexicana.combined.anc.bed',data.table = FALSE)
mex<-mexraw
colnames(mex)<-c('seqnames','start','end','testsnp','score')
mex$seqnames<-paste('chr',mex$seqnames,sep="")
#find the X% break
twop<-round(0.06*dim(mex)[1])

#increasing
cutoff<-mex$score[order(mex$score)[twop]]
#subset by score (mean mex ancestry)
mex<-mex[which(mex$score<cutoff),]

#decreasing
#cutoff<-mex$score[order(mex$score,decreasing = TRUE)[twop]]
#subset by score (mean mex ancestry)
#mex<-mex[which(mex$score>cutoff),]


#cM spline for maiz/mex

subsets<-lapply(1:10,function(cr){
  maiz_subset<-subset(maiz,seqnames==paste('chr',cr,sep=""))
  map_subset<-subset(map,chr==cr)
  maiz_subset$pos_cm <-spline(x=map_subset$pos_bp,y=map_subset$pos_cM,xout=maiz_subset$testsnp, method = "hyman")$y
  return(maiz_subset)
})

MaizNew<-bind_rows(subsets)
MaizNew$width<-MaizNew$end-MaizNew$start

grMaiz <- makeGRangesFromDataFrame(MaizNew, keep.extra.columns=TRUE)


#be aware that based on filtering you may not have hits on a chromosome and that will trip up things downstream
subsets<-lapply(c(1:10),function(cr){
  mex_subset<-subset(mex,seqnames==paste('chr',cr,sep=""))
  map_subset<-subset(map,chr==cr)
  mex_subset$pos_cm <-spline(x=map_subset$pos_bp,y=map_subset$pos_cM,xout=mex_subset$testsnp, method = "hyman")$y
  return(mex_subset)
})

MexNew<-bind_rows(subsets)
MexNew$width<-MexNew$end-MexNew$start


grMex <- makeGRangesFromDataFrame(MexNew, keep.extra.columns=TRUE)

#now compare overlaps

#Get an estimate of the sizes of tests that I'm looking at to get an expectation based on percentage


hitsMex <- countOverlaps(grRaisd, grMex, ignore.strand=TRUE)
#I don't think I want sum, I think I want is 0 or not
length(which(hitsMex>0))
raisdlen<-length(grRaisd)
raisdlenbases<-sum(data.frame(grRaisd)$width)
Mexlen<-length(grMex)
sum(data.frame( GenomicRanges::intersect(grRaisd, grMex,ignore.strand=TRUE))$width)

hitsMaiz <- countOverlaps(grRaisd, grMaiz, ignore.strand=TRUE)
length(which(hitsMaiz>0))
maizlen<-length(grMaiz)
sum(data.frame( GenomicRanges::intersect(grRaisd, grMaiz,ignore.strand=TRUE))$width)

hitsMexMaiz <- countOverlaps(grMex, grMaiz, ignore.strand=TRUE)
length(which(hitsMexMaiz>0))
sum(data.frame( GenomicRanges::intersect(grMex, grMaiz,ignore.strand=TRUE))$width)

#get approx proportions of bins from final:
proportions<-table(raisd$bin)/nrow(raisd)

bin<-as.numeric(cut_number(MaizNew$pos_cm,10))
MaizNew$bin<-bin
maiz<-MaizNew
maizproportions<-table(maiz$bin)

bin<-as.numeric(cut_number(MexNew$pos_cm,10))
MexNew$bin<-bin
mex<-MexNew
mexproportions<-table(mex$bin)

#Right now I'm not doing anything with rasid/mex or mex/maize overlaps but I'm calculating them just in case
#Not certain what the expectation or biological explanation would be if there were enrichment

#now iterate the function to sample randomly and build a bootstrap distribution

#also update raisdraw seqnames once to save time:
#raisdraw$seqnames<-paste('chr',raisdraw$seqnames,sep="")

#get the window sizes of raisd:
rdf<-data.frame(grRaisd)


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



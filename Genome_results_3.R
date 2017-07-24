#plink subamples

setwd("C:/Users/froug/Desktop/Chapter 3/Plink Subamples/Subample Results")

#subet only .genome files

b<-list.files(pattern="\\.genome")

#b<-b[grep("832", b)]

q<-list(0)

for (i in 1:length(b)){
  
  genome<-b[[i]]
  x<-read.csv(genome, header=TRUE, sep="", stringsAsFactors=FALSE)
  q[[i]]<-x
  names(q)[i]<-genome  
  #trimmed<-x[,10]
  #names(trimmed)<-genome
  #start<-cbind(start,trimmed)
}

#wait how do we deal with different numbers of dolphins!!!
#might have to always order them alphabetically

stset<-get("100snps275n1.genome",q)
start<-stset[,c(2,4)]
start <- apply(start,1,sort,decreasing=F)
start<-as.data.frame(t(start))
names(start)<-c("IID1", "IID2")

for (j in 1:length(q)){
  results<-q[[j]]
  
  for (k in 1:nrow(results)){
  if(results$IID1[k]>results$IID2[k]){
    holder<-results$IID1[k]
    results$IID1[k]<-results$IID2[k]
    results$IID2[k]<-holder
  }
  }
  
  start<-merge(start, results, by=c("IID1", "IID2"), all=TRUE)
  yy<-grep("pi_hat", ignore.case=TRUE, names(start))
  
  start<-start[,c(1:2, yy)]
  
  names(start)[ncol(start)]<-paste0(b[j],"_pi_hat")
  
}

#how to sort
View(start[order(start$PI_HAT, decreasing=TRUE),])

write.csv(start, "full_genome_pihat_results_final.csv")

#repeat for DST

stset<-get("100snps275n1.genome",q)
start<-stset[,c(2,4)]
start <- apply(start,1,sort,decreasing=F)
start<-as.data.frame(t(start))
names(start)<-c("IID1", "IID2")

for (j in 1:length(q)){
  results<-q[[j]]
  
  for (k in 1:nrow(results)){
    if(results$IID1[k]>results$IID2[k]){
      holder<-results$IID1[k]
      results$IID1[k]<-results$IID2[k]
      results$IID2[k]<-holder
    }
  }
  
  start<-merge(start, results, by=c("IID1", "IID2"), all=TRUE)
  yy<-grep("dst", ignore.case=TRUE, names(start))
  
  start<-start[,c(1:2, yy)]
  
  names(start)[ncol(start)]<-paste0(b[j],"_dst")
  
}

write.csv(start, "full_genome_dst_results_final.csv")


#plotting

IBD<-read.csv("full_genome_pihat_results_final.csv", row.names = 1)
IBS<-read.csv("full_genome_dst_results_final.csv", row.names = 1)

library(reshape2)

IBDm<-melt(IBD, id.vars=c("IID1", "IID2"))
IBDm$variable<-as.character(IBDm$variable)

x<-strsplit(IBDm$variable,c("X"),fixed=TRUE)
IBDm$dataset<-sapply(x, "[", 2)

x<-strsplit(IBDm$dataset,c("snps"),fixed=TRUE)
IBDm$snps<-sapply(x, "[", 1)

x<-strsplit(sapply(x, "[", 2),c("n"),fixed=TRUE)
IBDm$nind<-sapply(x, "[", 1)

#don't think you need iteration number but could pull that too

#add true values to data frame

safecopyIBDm<-IBDm

truth<-subset(IBDm, snps==nsnps & nind==nmax)[,c(1:2,4,6:7)]

truth2<-truth[!duplicated(truth),]

#truth2 is lookup table now do merge

IBDm2<-merge(IBDm, truth2, by=c("IID1", "IID2"), all.x=TRUE)

write.csv(IBDm2, "IBDm2.csv")

#remove NAs from value.x

IBDcom<-IBDm2[complete.cases(IBDm2),]

IBDcom$diff<-abs(IBDcom$value.y-IBDcom$value.x)

IBDcom$percenterr #what about when true value is 0? 

#average diffs across cat
#try summarize or aggregate?

testw<-aggregate(diff~snps.x+nind.x, data=IBDcom, mean)

safecopy<-testw
x<-safecopy

#then plot as a heatmap

x$snps.x<-as.factor(x$snps.x)
x$nind.x<-as.factor(x$nind.x)

p <-  ggplot(x, aes(x$snps.x, x$nind.x)) + 
  geom_tile(aes(fill = diff),colour = "white") +
  scale_fill_gradient(low = "green", high = "red", name="Difference") +
  labs(x="Number of SNPs", y="Number of Individuals")+
  ggtitle("Relatedness (IBD) Values")


p

##Repeat with IBS!!!

IBS<-read.csv("full_genome_dst_results_final.csv", row.names = 1)

library(reshape2)

IBSm<-melt(IBS, id.vars=c("IID1", "IID2"))
IBSm$variable<-as.character(IBSm$variable)

x<-strsplit(IBSm$variable,c("X"),fixed=TRUE)
IBSm$dataset<-sapply(x, "[", 2)

x<-strsplit(IBSm$dataset,c("snps"),fixed=TRUE)
IBSm$snps<-sapply(x, "[", 1)

x<-strsplit(sapply(x, "[", 2),c("n"),fixed=TRUE)
IBSm$nind<-sapply(x, "[", 1)

#don't think you need iteration number but could pull that too

#add true values to data frame

safecopyIBSm<-IBSm

nsnps=3863
nmax=275

truth<-subset(IBSm, snps==nsnps & nind==nmax)[,c(1:2,4,6:7)]

truth2<-truth[!duplicated(truth),]

#truth2 is lookup table now do merge

IBSm2<-merge(IBSm, truth2, by=c("IID1", "IID2"), all.x=TRUE)

write.csv(IBSm2, "IBSm2.csv")

#remove NAs from value.x

IBScom<-IBSm2[complete.cases(IBSm2),]

IBScom$diff<-abs(IBScom$value.y-IBScom$value.x)

IBScom$percenterr #what about when true value is 0? 

#average diffs across cat
#try summarize or aggregate?

testw<-aggregate(diff~snps.x+nind.x, data=IBScom, mean)

safecopy<-testw
x<-safecopy

#then plot as a heatmap

x$snps.x<-as.factor(as.numeric(x$snps.x))
x$nind.x<-as.factor(as.numeric(x$nind.x))

p <-  ggplot(x, aes(x$snps.x, x$nind.x)) + 
  geom_tile(aes(fill = diff),colour = "white") +
  scale_fill_gradient(low = "green", high = "red", name="Difference") +
  labs(x="Number of SNPs", y="Number of Individuals")+
  ggtitle("Genomic Similarity (IBS) Values")


p

#error distributions, overlayed
#make another one to show the variance?
#show max difference as well


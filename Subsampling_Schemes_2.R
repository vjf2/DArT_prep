setwd("C:/Users/froug/Desktop/DArT Seq/2017 results")

options(stringsAsFactors = FALSE)

ID_key<-read.csv("C:/Users/froug/Desktop/DArT Seq/2017 results/dolphin_sample_key.csv", colClasses = "character")

plink<-read.csv("plinked_filtered_2017.csv")

coancestry<-read.csv("franzinput_CR95.csv")

#set up subsampling scheme

#random draws individuals and snps

nsnps<-(ncol(coancestry)-1)/2
nmax<-nrow(coancestry)  

snpset<-c(50,100,200,400,800,1600,3200,nsnps)
nind<-c(10,20,50,100,200, nmax)

#####PLINK

plink<-plink[,(24:ncol(plink))] 

plink<-t(plink)

x<-merge(plink, ID_key[,1:2], by.x="row.names", by.y="Sample_ID")

x[,1]<-x$Dolphin_ID
pl<-x[,-ncol(x)]

#set up ped and map file

pl1<-data.frame(dummy_group=rep(1,nrow(pl)),pl)

m<-matrix(rep(0, 4*nrow(pl)), nrow=nrow(pl))

pl2<-data.frame(pl1[,1:2],m, pl1[,3:ncol(pl1)]) 

original.ped<-pl2

# have to add allele names to ped file 

a1<-paste0(1:nsnps, ".1")
a2<-paste0(1:nsnps, ".2")
aa<-sort(as.numeric(c(a1,a2)))
ag<-paste0(rep("G", nsnps), aa)

names(original.ped)[7:ncol(original.ped)]<-ag
  
original.map<-data.frame(rep(0,nsnps), paste0(rep("G", nsnps), 1:nsnps), rep(0, nsnps), rep(0, nsnps))

write.table(original.ped, "original.ped", sep="\t", col.names=FALSE, row.names=FALSE,quote=FALSE)
write.table(original.map, "original.map", sep="\t", col.names=FALSE, row.names=FALSE,quote=FALSE)

#source function at bottom to sample snps

nreps<-3

combos<-expand.grid(list(snpset,nind))

#da megaloop

type<-original.ped #or coancestry

#sample nind then put in snp_sample

for (i in 1:nrow(combos)){
  
  ns<-combos[i,1]
  inds<-combos[i,2]
  
  snp_sample(ns,inds,nreps,original.map, original.ped)

  }

#looks great! get a batch script for plink results



#########COANCESTRY



#let's do the coancestry subsamples first since those need to be done one-by-one

coancestry[1:10,1:10]

#read in dolphin sample key

x<-merge(coancestry, ID_key[,1:2], by.x="X", by.y="Sample_ID")

x[,1]<-x$Dolphin_ID
ca<-x[,-ncol(x)]

#draw random numbers of individuals

setwd("C:/Users/froug/Desktop/Chapter 3/Coancestry Subsamples")

#how many reps of each do you want? 

nreps<-3

combos<-expand.grid(list(snpset,nind))

#da megaloop

type<-ca #or plink

for (i in 1:nrow(combos)){
  
  nsnps<-combos[i,1]
  inds<-combos[i,2]
  
  ca_ind<-ca[sample(nrow(ca),inds),]

  ca_ind<-ca_ind[-1,]
  write.csv(ca_ind, file=paste0("S",nsnps,"I",inds,".csv"), row.names = FALSE)
  
  }

#ok now label the alleles so they can be subsetted properly


#plink function

snp_sample<-function(num_snps, num_inds, num_samp,original_map, original_ped){
  
  inds<-num_inds
  
  maps<-list(0)
  ped<-list(0)
  
  for (i in 1:num_samp) {
    x<-sample(original_map[,2], num_snps)
    new_map<-subset(original_map, original_map[,2] %in% x)
    
    snps1<-paste0(x,".1")
    snps2<-paste0(x,".2")
    snpnames<-c(snps1,snps2)
    prefix<-original_ped[,1:6]
    y<-colnames(original_ped)
    sampsnps<-original_ped[, which(colnames(original_ped) %in% snpnames)]
    new_ped<-cbind(prefix, sampsnps)
    
    new_ped<-new_ped[sample(nrow(new_ped),inds),]
    
    map_name=paste0(num_snps,"snps", inds, "n", i,".map")
    ped_name=paste0(num_snps,"snps", inds, "n", i,".ped")
    
    write.table(new_map, map_name, sep="\t", col.names=FALSE, row.names=FALSE,quote=FALSE)
    write.table(new_ped, ped_name, sep="\t", col.names=FALSE, row.names=FALSE,quote=FALSE)
    
  }
  
  
}


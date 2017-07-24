#V1
#Vivienne Foroughirad
#June 13, 2017 Created

#Read in raw data from SNPs and filter based on MAF, CallRate, Errors, Distance
#Format for Cervus, Plink, Coancestry, and Franz

#Read in raw datasets for single and double file, and ID key
setwd()

options(stringsAsFactors = FALSE)
library(HardyWeinberg)

single_raw<-read.csv(dart_output_singlerow, colClasses = "character", skip=6)
double_raw<-read.csv(dart_output_doublerow, colClasses = "character", skip=6)

#ID_key contains "Sample_ID", "Animal_ID", and optionally "Birthyear" and "Sex"

ID_key<-read.csv(ID_key, colClasses = "character")

mincallrate<-0.95 #default

#Check percent identity between duplicate samples, and remove snps from list that typed differently

#Identify duplicate samples


#dups is character vector with Animal_ID of duplicate sample
#assumes 5 duplicate samples, modify indices for more or less

dups<-c("ABC", "DEF", "GHI", "JKL", "MNO")
error_rate<-list()

#This calculates the error rates in duplicate samples

for (i in 1:length(dups)){

cdup<-ID_key$Sample_ID[ID_key$Animal_ID==dups[i]]
ctitle<-paste0("remove",i)

single_raw[,ctitle]<-ifelse(single_raw[,cdup[1]]==single_raw[,cdup[2]], "agree", "disagree")
error_rate[[i]]<-length(single_raw[,ctitle][single_raw[,ctitle]=="disagree"])/dim(single_raw)[1]
}

ers<-unlist(error_rate)

mean(ers)

#Sample error, remove duplicates

cdups<-ID_key$Sample_ID[ID_key$Animal_ID %in% dups]

minus<-cdups[6:10]

single<-single_raw[,!names(single_raw) %in% minus]
double<-double_raw[,!names(double_raw) %in% minus]

#actually remove the disagreeing ones

single<-subset(single,  single$remove1=="agree" &
                        single$remove2=="agree" &
                        single$remove3=="agree" &
                        single$remove4=="agree" &
                        single$remove5=="agree"
                 )

single<-single[,-((dim(single)[2]-4):(dim(single)[2]))]

#Calculate HWE for all alleles
#Convert reference/snp allele to major/minor allele

single$MN<-apply(single[,22:296], 1, function(x) length(x[x==2]))
single$NN<-apply(single[,22:296], 1, function(x) length(x[x==1]))
single$MM<-apply(single[,22:296], 1, function(x) length(x[x==0]))

single$maxA<-apply(single[,c("NN", "MM")], 1, max)

single$total<-apply(single[,c("MN","NN", "MM")], 1, sum)

single$MAF<-1-(single[,"maxA"]+((single[,"MN"])/2))/(single[,"total"])

single$MAcount<-single$total-single$maxA-single$MN

#Loop through and pass one value at a time to AA, AB, BB

single$pval<-rep(NA, dim(single)[1])
single$expectedAA<-rep(NA, dim(single)[1])
single$expectedAB<-rep(NA, dim(single)[1])
single$expectedBB<-rep(NA, dim(single)[1])

#Will give you warnings for expected counts below 5, but
#we're filtering out those anyway

for (i in 1:dim(single)[1]){
  AA<-single$maxA[i]
  AB<-single$MN[i]
  BB<-single$MAcount[i]
  
  output<-HWChisq(c(AA=AA, AB=AB, BB=BB), verbose=FALSE)
  
  single$pval[i]<-output$pval
  single$expectedAA[i]<-output$expected[1]
  single$expectedAB[i]<-output$expected[2]
  single$expectedBB[i]<-output$expected[3]
  
}

#Remove HWE, MAF, and CallRate cutoffs

single<-subset(single,  pval>=0.05 &
                        MAF>=0.01 &
                        CallRate>=mincallrate)

#Pull out one SNP per contig

single<-subset(single, Chrom_Tursiops_v14!="")

#Remove duplicates based on condition

#pick unique value from Chrom_Tursiops_v14 based on highest value of MAF

single<-single[with(single, ave(MAF, Chrom_Tursiops_v14, FUN=max)==MAF),]

#remove ties with equal highest MAF

single<-single[!duplicated(single$Chrom_Tursiops_v14),]

#now that we have the final set for single, match it up to double

#format Allele ID column for matching by creating unique ID

AlleleNumber<-strsplit(single$AlleleID,"\\|")
AlleleNumber<-unlist(lapply(AlleleNumber, "[[",1))

single<-cbind(single, AlleleNumber)

AlleleNumber2<-strsplit(double$AlleleID,"\\|")
AlleleNumber2<-unlist(lapply(AlleleNumber2, "[[",1))

double<-cbind(double, AlleleNumber2)

#filter double

double<-double[which(double$AlleleNumber2 %in% single$AlleleNumber),]

final_double<-subset(double, double$AlleleNumber2 %in% single$AlleleNumber)

#still duplicates because of multiple snps in the same read
#give double a dummy ID

n<-dim(double)[1]/2

double$allele_unique<-rep(1:n, each=2)

final_ids<-subset(double$allele_unique, double$AlleleID %in% single$AlleleID)

final_double<-subset(double, double$allele_unique %in% final_ids)

#write.csv(final_double, "check_double.csv")

#Reformat final double for CERVUS, plink, FRANZ, etc. 

cervus_double<-t(final_double)

n2<-dim(cervus_double)[2]/2

al<-rep(c("a","b"), n2)

cervus_double<-as.data.frame(rbind(cervus_double, al))

tt<-paste0(cervus_double["AlleleNumber2",], cervus_double["al",])

names(cervus_double)<-tt

#write.csv(cervus_double, "cervus_double_noMAF.csv")
#post processing
#convert - to * and 0 to 2
#remove excess rows

#Format the data for Plink, Franz
#plink is letter, franz is 1-4 based on letter

plink_double<-t(cervus_double)

#Make nuc1 and nuc0 columns

lt<-as.data.frame(plink_double[,c("SNP","allele_unique")], na.strings=c("", NA))

names(lt)<-c("SNP1", "allele_unique")

lt$SNP1[lt$SNP1==""] <- NA

lt<-lt[complete.cases(lt),]

plink_double<-merge(plink_double, lt, by="allele_unique")

AlleleLetter<-strsplit(plink_double$SNP1,":")
plink_double$AlleleLetter<-unlist(lapply(AlleleLetter, "[[",2))

RefLetter<-strsplit(plink_double$AlleleLetter, ">")
plink_double$RefLetter<-unlist(lapply(RefLetter, "[[",1))
plink_double$AltLetter<-unlist(lapply(RefLetter, "[[",2))

plink_double$Nuc0<-ifelse(plink_double$al=="b", plink_double$AltLetter, plink_double$RefLetter)
plink_double$Nuc1<-ifelse(plink_double$al=="a", plink_double$AltLetter, plink_double$RefLetter)

x<-plink_double

for (i in 23:297) {
  
  x[x[,i]=="0", names(x)[i]]<-x[x[,i]=="0", "Nuc0"]
  x[x[,i]=="1", names(x)[i]]<-x[x[,i]=="1", "Nuc1"]
  
}

x[,23:297]<-apply(x[,23:297], 2, function(x) gsub("-", 0, x))

plink_filtered_2017<-x

#write.csv(plink_filtered_2017, "plinked_filtered_2017.csv")

#Franz and coancestry convert letters to numbers

x<-plink_filtered_2017

for (i in 23:297) {
  
  #lt$SNP1[lt$SNP1==""] <- NA
  
  x[,i][x[,i]=="A"]<-1
  x[,i][x[,i]=="T"]<-2
  x[,i][x[,i]=="G"]<-3
  x[,i][x[,i]=="C"]<-4
  
}

x[,23:297]<-apply(x[,23:297], 2, function(x) gsub("-", 0, x))

franz<-x[,(23:297)]
franz<-t(franz)

write.csv(franz, "franzinput_CR95.csv")


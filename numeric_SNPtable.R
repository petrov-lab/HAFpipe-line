args<-commandArgs(TRUE)

### libraries
library(data.table)

## FILES
snpFile <- args[1] #snpFile<-"/mnt/cages/ref_data/dgrp/snps/all_lines/new_2L.csv"
#snpInfoFile <- args[2] #snpInfoFile<-"/mnt/cages/ref_data/dgrp/snps/all_lines/new_2L.csv.altCt" 
snpInfoFile <- paste0(snpFile,".alleleCts")
###### snpInfoFile must be in the same row-order as snpFile, and have a column for identity of the alt allele

## READ IN DATA
SNPs<- as.data.frame(fread(snpFile, header=TRUE))
snpInfo<- as.data.frame(fread(snpInfoFile))
chrom<-colnames(SNPs)[1]

##MAKE NUMERIC SNP TABLE
numericSNPs<-apply(SNPs[,3:dim(SNPs)[2]],2,function(x) .5*(as.numeric(x=="N")) + as.numeric(x==snpInfo$alt))
rownames(numericSNPs)<-SNPs[,1]
headers<-paste(c(chrom,colnames(numericSNPs)),collapse=",")
cat(paste(headers,'\n',sep=""),file=paste(snpFile,".numeric",sep=""))
write.table(numericSNPs,file=paste(snpFile,".numeric",sep=""),append=TRUE,quote=FALSE,sep=",",row.names=TRUE,col.names=FALSE)


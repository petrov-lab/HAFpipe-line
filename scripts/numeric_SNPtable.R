args<-commandArgs(TRUE)

# DEPRECATED
# This script reads the whole tables into memory at once, which is a problem in itself,
# but also for larger datasets fails even with enough memory (~150GB in our tests),
# as we reach the limit of 2^31-1 values that R can handle in a single table...
# See https://github.com/petrov-lab/HAFpipe-line/issues/7 for details.
# We hence re-implemented this step in Python - use the .py script instead!

### libraries
library(data.table)

## FILES
snpFile <- args[1]
snpInfoFile <- paste0(snpFile,".alleleCts")
###### snpInfoFile must be in the same row-order as snpFile, and have a column for identity of the alt allele

## READ IN DATA
SNPs<- as.data.frame(fread(snpFile, header=TRUE))
snpInfo<- as.data.frame(fread(snpInfoFile))
chrom<-colnames(SNPs)[1]

##MAKE NUMERIC SNP TABLE
### CODES: ref=0; alt=1; het=.5; missing=-1

make_numeric=function(x){-0.5*(as.numeric(x==snpInfo$ref)) + 0.5*as.numeric(x==snpInfo$alt) + -1.5*as.numeric(x=="N") + .5}
numericSNPs<-apply(SNPs[,3:dim(SNPs)[2]],2,make_numeric)

### OLD: numericSNPs<-apply(SNPs[,3:dim(SNPs)[2]],2,function(x) {.5*(as.numeric(x=="N")) + as.numeric(x==snpInfo$alt))
rownames(numericSNPs)<-SNPs[,1]
headers<-paste(c(chrom,colnames(numericSNPs)),collapse=",")
cat(paste(headers,'\n',sep=""),file=paste(snpFile,".numeric",sep=""))
write.table(numericSNPs,file=paste(snpFile,".numeric",sep=""),append=TRUE,quote=FALSE,sep=",",row.names=TRUE,col.names=FALSE)

args<-commandArgs(TRUE)

### libraries
library(data.table)
library(foreach)

## READ
snpFile <- args[1]
pctMissing <- as.numeric(args[2])
snps=fread(snpFile)

##ADD Ns
nGenotypes=(ncol(snps)-2)*nrow(snps)
nMissing=nGenotypes*pctMissing/100
ixMissing=cbind(sample(1:nrow(snps),nMissing,replace=TRUE),sample(3:ncol(snps),nMissing,replace=TRUE))
nn=foreach(ii=3:ncol(snps))%do%{mask=ixMissing[,2]==ii;sites=unique(ixMissing[mask,1]);snps[sites,ii]="N";ii}

## WRITE
write.csv(snps,file=paste0(snpFile,".",pctMissing,"pctMissing"),quote=F,row.names=F)
cat(round(1000*sum(snps=="N")/nGenotypes)/10, "% missing calls added and written to ",paste0(snpFile,".",pctMissing,"pctMissing"))
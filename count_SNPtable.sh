#!/bin/bash

## calculate allele frequency for founders at each site in SNP table  #
########################
snpTable=${1}
#########################
echo "pos,ref,alt,nCt,refCt,altCt,hetCt" > ${snpTable}.alleleCts
cat $snpTable | awk -F ',' '
(NR>1){ 
	pos=$1; ref=$2; alt="o"
	nCt=gsub("N","N",$0)
	refCt=gsub(ref,ref,$0)-1
	altCt=0;hetCt=0;
	
	for (i=3;i<=NF;i++) {
		if($i!="N" && $i!=ref){
			if(match("ACGT",$i)>0){ alt=$i; altCt++; }
			else {hetCt++}
		}
	}
	if ((altCt+nCt+refCt+hetCt+2)!=NF) { print "invalid counts! "pos","ref","alt","nCt","refCt","altCt","hetCt; exit 1 } 
	print pos","ref","alt","nCt","refCt","altCt","hetCt 
}' >> ${snpTable}.alleleCts

echo "allele counts written to:"
echo ${snpTable}.alleleCts
			

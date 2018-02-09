#!/bin/bash

## calculate allele frequency for founders at each site in SNP table  #
########################
#snpTable=/mnt/cages/ref_data/dgrp/snps/all_lines/noN/new_2L.csv
snpTable=${1}
#########################
echo "allele freqs will be written to ${snpTable}.alleleCts"
echo "pos,ref,alt,nCt,refCt,altCt" > ${snpTable}.alleleCts
cat $snpTable | awk -F ',' '
(NR>1){ 
	pos=$1; ref=$2; alt="o"
	nCt=gsub("N","N",$0)
	refCt=gsub(ref,ref,$0)-1
	altCt=0;
	#evenAF=0; propAF=0; nCt=0;
	for (i=3;i<=NF;i++) {
		if($i!="N" && $i!=ref)
		{ alt=$i; altCt++; }
	#	weight=1;
	#	if ($i=="N") { weight=.5; nCt++ }
	#	else if ($i==ref) { weight=0 }
	#	else {alt=$i}
	#	evenAF+=weight
	#	propAF+=weight*fProp[i]
	}
	if ((altCt+nCt+refCt+2)!=NF) { print "invalid counts! "pos","ref","alt","nCt","refCt","altCt; exit 1 } 
	print pos","ref","alt","nCt","refCt","altCt 
}' >> ${snpTable}.alleleCts

			

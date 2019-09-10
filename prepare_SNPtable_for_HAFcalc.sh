#!/bin/bash

## prepare SNP table for allele frequency calculation  #
########################
snptable=${1}
scriptdir=$(dirname "$0")
#########################

if [ ! -e ${snptable}.numeric ]; then
	echo "making numeric version of $snptable"; 
	Rscript $scriptdir/numeric_SNPtable.R $snptable
fi				

echo "bgzipping ${snptable}.numeric"

chrom=$(head -1 $snptable | cut -f1 -d',')
tail -n +2 ${snptable}.numeric | tr ',' '\t' | awk -v chrom="$chrom" '{
		alleleCt=NF-1
		printf("%s\t%s",chrom,$1 ); $1="";

		calledAlts=gsub(1,1,$0); uncalled=gsub(5,5,$0); 
		calledRate=calledAlts/(alleleCt-uncalled);  
		gsub("0.5",calledRate,$0); 
		print
	}' | bgzip > ${snptable}.numeric.bgz
tabix -s 1 -b 2 -e 2 ${snptable}.numeric.bgz


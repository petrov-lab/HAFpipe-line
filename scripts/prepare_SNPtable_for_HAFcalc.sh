#!/bin/bash

## prepare SNP table for allele frequency calculation by
## assigning missing calls a fractional count equal to
## the frequency of alt alleles among called genotypes #
########################
snptable=${1}
maindir=$(dirname "$0")/..
#########################

if [ ! -e ${snptable}.numeric ]; then
    echo "making numeric version of $snptable";
    Rscript ${maindir}/scripts/numeric_SNPtable.R $snptable
fi

echo "bgzipping ${snptable}.numeric"

chrom=$(head -1 $snptable | cut -f1 -d',')
tail -n +2 ${snptable}.numeric | tr ',' '\t' | awk -v chrom="$chrom" '{
        alleleCt=NF-1
        printf("%s\t%s",chrom,$1 );
        $1="";
        altCt=0; calledCt=0
        for(ii=2;ii<=NF;ii++){
            if($ii>=0){
                altCt=altCt+$ii;
                calledCt=calledCt+1
            }
        }
        calledRate=altCt/calledCt;
        for(ii=2;ii<=NF;ii++){
            if($ii<0){
                $ii=calledRate
            }
        }
        print
    }' | bgzip > ${snptable}.numeric.bgz
tabix -s 1 -b 2 -e 2 ${snptable}.numeric.bgz

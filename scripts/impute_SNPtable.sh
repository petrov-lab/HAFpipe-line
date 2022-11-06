#!/bin/bash

# ==================================================================================================
#      IMPUTE ALLELES IN SNP FILE
# ==================================================================================================

maindir=$(dirname "$0")/..

snpFile=${1}
nofFounders=$(( $(head -1 $snpFile | tr ',' ' ' | wc -w) - 2 ))

if [ ! -f ${snps}.alleleCts ]; then
    echo "counting alleles in $snps"
    ${maindir}/scripts/count_SNPtable.sh $snps
fi

posCol=$(head -1 $snpFile.alleleCts | tr ',' '\n' | awk '($1=="pos"){print NR}')
refCtCol=$(head -1 $snpFile.alleleCts | tr ',' '\n' | awk '($1=="refCt"){print NR}')
altCtCol=$(head -1 $snpFile.alleleCts | tr ',' '\n' | awk '($1=="altCt"){print NR}')
altCol=$(head -1 $snpFile.alleleCts | tr ',' '\n' | awk '($1=="alt"){print NR}')
namesRow=$(head $snpFile | grep -n Ref | cut -f1 -d':')

head -$namesRow $snpFile | tail -1 > $snpFile.imputed

join -j 1 \
<(tail -n +2 $snpFile.alleleCts | awk -v cols="$posCol,$refCtCol,$altCtCol,$altCol" -F ',' '
    BEGIN{split(cols,colArr,",")}{print $colArr[1]"\t"$colArr[2]"\t"$colArr[3]"\t"$colArr[4]}' ) \
<(tail -n +$(( 1 + $namesRow )) $snpFile  | tr ',' '\t' ) | \
awk -F ' ' -v nofFounders="$nofFounders" '
BEGIN{ srand()}
{
    pos=$1; refCt=$2; altCt=$3; alt=$4; ref=$5
    if ((refCt+altCt)>0){
    printf pos","ref
    for (i=6;i<=NF;i++) {
        if ($i=="N")
        {   if ( rand() <= altCt/(refCt+altCt))
                { printf ","alt }
            else {printf ","ref }
        }
        else {printf ","$i}
    }
    printf"\n"
    }
}'  >> $snpFile.imputed
echo "imputed SNPs written to
${snpFile}.imputed"

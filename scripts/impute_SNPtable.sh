#!/bin/bash

# ==================================================================================================
#      IMPUTE ALLELES IN SNP FILE
# ==================================================================================================

maindir=$(dirname "$0")/..

snptable=${1}
nofFounders=$(( $(head -1 ${snptable} | tr ',' ' ' | wc -w) - 2 ))

if [ ! -f ${snps}.alleleCts ]; then
    echo "counting alleles in $snps"
    ${maindir}/scripts/count_SNPtable.sh $snps
fi

posCol=$(head -1 ${snptable}.alleleCts | tr ',' '\n' | awk '($1=="pos"){print NR}')
refCtCol=$(head -1 ${snptable}.alleleCts | tr ',' '\n' | awk '($1=="refCt"){print NR}')
altCtCol=$(head -1 ${snptable}.alleleCts | tr ',' '\n' | awk '($1=="altCt"){print NR}')
altCol=$(head -1 ${snptable}.alleleCts | tr ',' '\n' | awk '($1=="alt"){print NR}')
namesRow=$(head ${snptable} | grep -n Ref | cut -f1 -d':')

head -$namesRow ${snptable} | tail -1 > ${snptable}.simpute

join -j 1 \
<(tail -n +2 ${snptable}.alleleCts | awk -v cols="$posCol,$refCtCol,$altCtCol,$altCol" -F ',' '
    BEGIN{split(cols,colArr,",")}{print $colArr[1]"\t"$colArr[2]"\t"$colArr[3]"\t"$colArr[4]}' ) \
<(tail -n +$(( 1 + $namesRow )) ${snptable}  | tr ',' '\t' ) | \
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
}'  >> ${snptable}.simpute
echo "imputed SNPs written to ${snptable}.simpute"

# We want to create the numeric table now, so that this is done before we move
# on to Tasks 3 and 4. If running these tasks in parallel for many samples,
# it could otherwise happen that we work with corrupted files, see
# https://github.com/petrov-lab/HAFpipe-line/issues/5
if [ ! -e ${snptable}.simpute.numeric.bgz ]; then
    if [ ! -e ${snptable}.simpute.alleleCts ]; then
        echo "counting alleles in ${snptable}.simpute"
        ${maindir}/scripts/count_SNPtable.sh ${snptable}.simpute
    fi

    echo "preparing ${snptable}.simpute for allele frequency calculation"
    ${maindir}/scripts/prepare_SNPtable_for_HAFcalc.sh ${snptable}.simpute
    if [ ! -e ${snptable}.simpute.numeric.bgz ]; then
        echo "creating ${snptable}.simpute.numeric failed"
        exit 1
    fi
fi

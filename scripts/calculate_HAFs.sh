#!/bin/bash

# ==================================================================================================
#      Usage
# ==================================================================================================

usage()
{
    echo "usage: calculate_HAFs.sh
                      [ -f freqfile ]
                      [ -s snptable ]
                      [ -o outdir=<dirname freqfile> ]
                      [ -d maindir=<main directory of HAF-pipe> ]
                      [ -h help ]
"
}
if [ $# -lt 1 ]; then usage; exit; fi

# ==================================================================================================
#      Command Line Arguments
# ==================================================================================================

maindir=$(dirname "$0")/..

while [ "$1" != "" ]; do
    case $1 in
        -f | --freqfile )       shift
                                freqs=$1
                                ;;
        -s | --snptable )       shift
                                snptable=$1
                                ;;
        -o | --outdir )         shift
                                outdir=$1
                                ;;
         -d | --maindir )       shift
                                maindir=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     echo unknown flag $1 ; usage
                                exit 1
    esac
    shift
done

# ==================================================================================================
#      Main
# ==================================================================================================

outfile=$outdir/$(basename $freqs | sed 's/.freqs/.afSite/')
chrom=$(head -1 $freqs | cut -f1 -d' ')

echo "snps are: $snptable"
#ls -lh ${snptable}.bgz.tbi
#ls -lh ${snptable}.numeric
#ls -lh ${snptable}.alleleCts

if [ ! -e ${snptable}.numeric.bgz ]; then
    if [ ! -e ${snptable}.alleleCts ]; then
        echo "counting alleles in $snptable"
        ${maindir}/scripts/count_SNPtable.sh $snptable
    fi
    echo "preparing $snptable for allele frequency calculation"
    ${maindir}/scripts/prepare_SNPtable_for_HAFcalc.sh $snptable
fi

echo "pos,af" > $outfile
cat $freqs | tr ' ' '\t' | awk -v snptable="${snptable}.numeric.bgz" '{
    print $0;
    system("tabix "snptable" "$1":"$2"-"$3" | cut -f2-")
}' | awk -v chrom="$chrom" '
    ($1==chrom){ for(ii=4;ii<=NF;ii++){freq[ii-3]=$ii} };
    ($1!=chrom){
        snp_win_AF=0
        for(ii=2;ii<=NF;ii++){
            snp_win_AF=snp_win_AF+$ii*freq[ii-1]
        }
        snp_sum_AF[$1]=snp_sum_AF[$1]+snp_win_AF
        snp_win_ct[$1]=snp_win_ct[$1]+1
    }
    END { for (pos in snp_sum_AF) print pos"\t"snp_sum_AF[pos]/snp_win_ct[pos] }
' | sort -k1n | tr '\t' ',' >>  $outfile

echo "HAFs written to $outfile"

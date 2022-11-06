#!/bin/bash

#calculate HAFs
usage()
{
    echo "usage: calculate_HAFs.floorHapFreqs.sh
                    [ -f freqfile ]
                    [ -s snptable ]
                    [ -p poolSize=100 ]
                    [ -o outdir=<dirname freqfile> ]
                    [ -d scriptdir=<dirname of this script> ]
                    [ -h help ]
"
}
if [ $# -lt 1 ]; then usage; exit; fi

# ==================================================================================================
#      Command Line Arguments
# ==================================================================================================

scriptdir=$(dirname "$0")
poolSize=100

while [ "$1" != "" ]; do
    case $1 in
        -f | --freqfile )       shift
                                freqs=$1
                                ;;
        -s | --snptable )       shift
                                snptable=$1
                                ;;
        -p | --poolSize )       shift
                                poolSize=$1
                                ;;
        -o | --outdir )         shift
                                outdir=$1
                                ;;
         -d | --scriptdir )      shift
                                scriptdir=$1
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

if [ -z $outdir ]; then outdir=$(dirname $freqs); fi
outfile=$outdir/$(basename $freqs | sed 's/.freqs/.afSite/')
chrom=$(head -1 $freqs | cut -f1 -d' ')

echo "snps are: $snptable"
#ls -lh ${snptable}.bgz.tbi
#ls -lh ${snptable}.numeric
#ls -lh ${snptable}.alleleCts

if [ ! -e ${snptable}.bgz ]; then
     if [ ! -e ${snptable}.numeric ]; then
        if [ ! -e ${snptable}.alleleCts ]; then
            echo "counting alleles in $snptable"; $scriptdir/count_SNPtable.sh $snptable
        fi
        echo "making numeric version of $snptable"; Rscript $scriptdir/numeric_SNPtable.R $snptable
    fi
    echo "bgzipping and indexing $snptable";
    tail -n +2 ${snptable}.numeric | tr ',' '\t' | awk -v chrom="$chrom" '{
        alleleCt=NF-1
        printf("%s\t%s",chrom,$1 ); $1="";

        calledAlts=gsub(1,1,$0); uncalled=gsub(5,5,$0);
        calledRate=calledAlts/(alleleCt-uncalled);
        gsub("0.5",calledRate,$0);
        print
    }' | bgzip > ${snptable}.bgz
    tabix -s 1 -b 2 -e 2 ${snptable}.bgz
fi
if [ ! -e ${snptable}.bgz.tbi ]; then tabix -s 1 -b 2 -e 2 ${snptable}.bgz; fi

echo "pos,af" > $outfile
cat $freqs | tr ' ' '\t' | awk -v snptable="${snptable}.bgz" '{
    print $0;
    system("tabix "snptable" "$1":"$2"-"$3" | cut -f2-")
}' | awk -v chrom="$chrom" -v poolSize="$poolSize" '
    ($1==chrom){
        totalfreq=0
        for(ii=4;ii<=NF;ii++){
            if($ii<(1/(poolSize*2))){$ii=0};
            if($ii>(1-(1/(poolSize*2)))){$ii=1};
            freq[ii-3]=$ii
            totalfreq=totalfreq+$ii
        }
    };
    ($1!=chrom){
        snp_win_AF=0
        for(ii=2;ii<=NF;ii++){
            snp_win_AF=snp_win_AF+($ii*freq[ii-1]/totalfreq)
        }
        snp_sum_AF[$1]=snp_sum_AF[$1]+snp_win_AF
        snp_win_ct[$1]=snp_win_ct[$1]+1
    }
    END { for (pos in snp_sum_AF) print pos"\t"snp_sum_AF[pos]/snp_win_ct[pos] }
' | sort -k1n | tr '\t' ',' >>  $outfile

echo "HAFs written to $outfile"

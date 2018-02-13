#!/bin/bash
	
#calculate HAFs
usage()
{
    echo "usage: calculate_HAFs.sh [ -f freqfile ] | [ -fl filename of freqfile list ]  
					  [ -s snps ] | [ -sl filename of snps list ]
					  [ -o outDir=<dirname freqfile> ] 
					  [ -t threads=1(default)] 
					  [ -h help ]
"
}
if [ $# -lt 1 ]; then usage; exit; fi

##### Main 
threads=1
scriptsDir=$HOME/scripts/HAFpipe

while [ "$1" != "" ]; do
    case $1 in
        -f | --freqfile )        shift
                                freqs=($1)
                                ;;
        -fl | --freqlist )       shift
                                freqs=($(cat $1))
                                ;;
        -s | --snps )           shift
                                snps=($1)
                                ;;
        -sl | --snpslist )      shift
                                snps=($(cat $1))
                                ;;
        -o | --outDir )         shift
                                outDir=$1
                                ;;
        -t | --threads )        shift
                                threads=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     echo unknown flag $1 ; usage
                                exit 1
    esac
    shift
done

chrom=$(head -1 ${snps}.numeric | cut -f1 -d',')


get_HAFs(){
	## VARS 
	snps=${1};
	outDir=${2}
	freqs=${3}

outFile=$outDir/$(basename $freqs | sed 's/.freqs/.afSite/')

if [ ! -f ${snps}.bgz.tbi ]; then
 	if [ ! -f ${snps}.numeric ]; then
		if [ ! -f ${snps}.alleleCts ]; then
			echo "counting alleles in $snps"; count_SNPtable.sh $snps	
		fi
		echo "making numeric version of $snps"; Rscript numeric_SNPtable.R $snps
	fi				
	echo "bgzipping and indexing $snps"; 
	tail -n +2 $snps.numeric | tr ',' '\t' | awk -v chrom="$chrom" '{
		printf("%s\t%s",chrom,$1 ); $1=""; 
		calledAlts=gsub(1,1,$0); uncalled=gsub(5,5,$0); calledRate=calledAlts/(NF-uncalled-1);  
		gsub("0.5",calledRate,$0); 
		print
	}' | bgzip > ${snps}.bgz
	tabix -s 1 -b 2 -e 2 ${snps}.bgz
fi

echo "pos,af" > $outFile
cat $freqs | tr ' ' '\t' | awk -v snps="${snps}.bgz" '{ 
	print $0; 
	system("tabix "snps" "$1":"$2"-"$3" | cut -f2-") 
}' | awk -v chrom="$chrom" '
	($1==chrom){ for(ii=4;ii<=NF;ii++){freq[ii]=$ii} }; 
	($1!=chrom){  
		snp_win_AF=0
		for(ii=2;ii<=NF;ii++){
			snp_win_AF=snp_win_AF+$ii*freq[ii-1]
		}
		snp_sum_AF[$1]=snp_sum_AF[$1]+snp_win_AF
		snp_win_ct[$1]=snp_win_ct[$1]+1
	}
	END { for (pos in snp_sum_AF) print pos"\t"snp_sum_AF[pos]/snp_win_ct[pos] }' | sort -k1n | tr '\t' ',' >>  $outFile

export -f get_HAFs
parallel --gnu -j${threads} get_HAFs $snps $outDir ::: ${freqs[*]}

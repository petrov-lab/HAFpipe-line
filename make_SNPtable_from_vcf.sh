#!/bin/bash

###############
## HELP
#################
usage()
{
    echo "usage: make_SNPtable_from_vcf.sh [-v --vcf] [-c --chrom] [ -o --outfile ] 
    optional: [-f --firstSampleCol (default 10)] [-m --mincalls (default 2)] [-k --keephets] [ -t --threads (default 1)] [ -h --help ]
	** note that the vcf may be gzipped 
	** only biallelic sites with at least one ref and one alt call will be recorded
	** use --mincalls to limit to sites with a higher number of called alleles
"
}
if [ $# -lt 1 ]; then usage; exit; fi

##### Main

firstSampleCol=10
threads=1
keephets=0

while [ "$1" != "" ]; do
    case $1 in
        -v | --vcf )            shift
                                vcf=($1)
                                ;;
        -c | --chrom )          shift
                                chrom=$1
                                ;;
        -f | --firstSampleCol ) shift
                                firstSampleCol=$1
                                ;;
        -m | --mincalls )       shift
                                mincalls=$1
                                ;;
        -k | --keephets )       
                                keephets=1
                                ;;                            
        -o | --outfile )        shift
                                outfile=$1
                                ;;
        -t | --threads )        shift
                                threads=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     echo "unknown flag $1"
				usage
                                exit 1
    esac
    shift
done

###################
## MAIN
##################

echo "making snptable for chrom $chrom from $vcf, starting from column $firstSampleCol, and writing to $outfile"

	
	##PRINT HEADER
	zcat -f $vcf | head -1000  | grep "#C" | head -1 | awk -v fsc=$firstSampleCol '{for(i=fsc;i<=NF;i++){printf $i","}}' |  sed "s/^/${chrom},Ref,/" | sed 's/,$/\n/'  > $outfile 
	
	zcat -f $vcf | grep -P '^'$chrom'\t' | grep PASS | awk -v fsc="$firstSampleCol" -v mincalls="$mincalls" -v keephets="$keephets" '
	BEGIN{
		baseCodes["AG"]="R";baseCodes["GA"]="R"
		baseCodes["CT"]="Y";baseCodes["TC"]="Y"
		baseCodes["CG"]="S";baseCodes["GC"]="S"
		baseCodes["AT"]="W";baseCodes["TA"]="W"
		baseCodes["GT"]="K";baseCodes["TG"]="K"
		baseCodes["AC"]="M";baseCodes["CA"]="M"
	}
	{
		if(substr($0, 0, 1)!="#") {
			if(length($5)==1 && length($4)==1) {
				refCt=gsub("0/0","0/0",$0)-1; 
				altCt=gsub("1/1","1/1",$0); 
				if(refCt>0 && altCt>0 && ((refCt+altCt) > minCalls)){
					
					printf $2","$4 
			
					for(i=fsc; i<=NF; i++) {
						split($(i),parts,":")
						GT = parts[1]
						printf "," 
						if(GT=="0/0") {
							printf $4 
						} else if (GT=="1/1") {
							printf $5 
						} else if (keephets>0 && ((GT=="0/1") || (GT=="1/0"))) {
							printf baseCodes[$4$5] 
						} 
						else {
							printf "N"
						} 
						
					}
					print "" 
				}
			}
		}		
	}' >> $outfile

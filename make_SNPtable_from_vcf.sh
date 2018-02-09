#!/bin/bash

###############
## HELP
#################
usage()
{
    echo "usage: make_SNPtable_from_vcf.sh [-v vcf.gz] [-c chroms] [ -o outDir ] [-f firstSampleCol=10] [-m minNofCalls=2] [ -t threads] [ -h help ]
	** note that the vcf must be gzipped **
"
}
if [ $# -lt 1 ]; then usage; exit; fi

##### Main

firstSampleCol=10
threads=1

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
        -m | --minCalls )       shift
                                minCalls=$1
                                ;;
        -o | --outdir )         shift
                                outDir=$1
                                ;;
        -t | --threads )        shift
                                threads=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

###################
## MAIN
##################
outfileBase=$(basename $vcf | sed 's/.vcf.gz//'); fi
outFile=${outDir}/${outfileBase}.${chrom}.snptable;

echo "making snptable for chroms: $chrom from $vcf, starting from column $firstSampleCol, and writing to ${outFile}"

	
	##PRINT HEADER
	zcat $vcf | head -1000  | grep "#C" | head -1 | awk -v fsc=$firstSampleCol '{for(i=fsc;i<=NF;i++){printf $i","}}' |  sed "s/^/${chrom},Ref,/" | sed 's/,$/\n/'  > $outFile 
	
	zcat $vcf | grep -P '^'$chrom'\t' | awk -v fsc="$firstSampleCol" -v minCalls="$minCalls" '
	{
		if(substr($0, 0, 1)!="#") {
			if(length($5)==1 && length($4)==1) {
				refCt=gsub("0/0","0/0",$0)-1; 
				altCt=gsub("1/1","1/1",$0); 
				if(refCt>0 && altCt>0 && (refCt+altCt) > minCalls){
					
					printf $2","$4 
			
					for(i=fsc; i<=NF; i++) {
						split($(i),parts,":")
						GT = parts[1]
						printf "," 
						if(GT=="0/0") {
							printf $4 
						} else if (GT=="1/1") {
							printf $5 
						} 
						else {
							printf "N"
						} 
						
					}
					print "" 
				}
			}
		}		
	}' >> $outFile

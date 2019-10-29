#!/bin/bash

###############
## HELP
#################
usage()
{
    echo "usage: make_SNPtable_from_vcf.sh [-v --vcf] [-c --chrom] 
    optional: 
    [ -s --snptable (default: $(echo $vcf | sed 's/.gz$//' | sed 's/.vcf//').$chrom.snptable)] 
    [-f --firstSampleCol (default 10) ] 
    [-m --mincalls (default 2) ] 
    [-u --subsetlist (default none)]
    [-k --keephets ] 
    [ -t --threads (default 1) ] 
    [ -h --help ]
	** note that the vcf may be gzipped 
	** only biallelic sites with at least one ref and one alt call will be recorded
	** use --mincalls to limit to sites with a higher number of called alleles
	** subset list is a 1-column list of names of founder haplotypes to keep in the snptable; 
	   if not supplied, all founders will be included
"
}
if [ $# -lt 1 ]; then usage; exit; fi

##### Main

firstSampleCol=10
threads=1
keephets=0
mincalls=2
subsetlist=none

while [ "$1" != "" ]; do
    case $1 in
        -v | --vcf )            shift
                                vcf=$1
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
        -u | --subsetlist )     shift
                                subsetlist=$1
                                ;;
        -s | --snptable )       shift
                                snptable=$1
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

if [ -z $snptable ]; then snptable=$(echo $vcf | sed 's/.gz$//' | sed 's/.vcf//').$chrom.snptable; fi

###################
## MAIN
##################

echo "making snptable for chrom $chrom from $vcf, starting from column $firstSampleCol"
scriptDir=$(dirname $0)
	
	##PRINT HEADER
	zcat -f $vcf | head -1000  | grep "#C" | head -1 | awk -v fsc=$firstSampleCol '{for(i=fsc;i<=NF;i++){printf $i","}}' |  sed "s/^/${chrom},Ref,/" | sed 's/,$/\n/'  > $snptable 
	
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
		if(length($5)==1 && length($4)==1) {
			refCt=gsub("0/0","0/0",$0); 
			altCt=gsub("1/1","1/1",$0); 
			hetCt=gsub("0/1","0/1",$0); 
			missingCt=gsub("\\./\\.","./.",$0); 
			
			if(refCt>0 && altCt>0 && ((refCt+altCt) > mincalls)){
				
				printf $2","$4 
		
				for(i=fsc; i<=NF; i++) {
					split($(i),parts,":")
					GT = parts[1]
					printf "," 
					if(GT=="0/0") {
						printf $4 
					} else if (GT=="1/1") {
						printf $5 
					} else if (keephets>0 && (GT=="0/1")) {
						printf baseCodes[$4$5] 
					} 
					else {
						printf "N"
					} 
					
				}
				print "" 
			}
		}		
	}' >> $snptable

	if [ ! "$subsetlist" == "none" ]; then
		$scriptDir/Extra/subset_SNPtable.sh \
		-s $snptable \
		-o ${snptable}.subset \
		-f $subsetlist \
		-m $mincalls ;
		mv ${snptable}.subset $snptable
	fi

	echo "SNP table written to:"
	echo $snptable

	$scriptDir/count_SNPtable.sh $snptable

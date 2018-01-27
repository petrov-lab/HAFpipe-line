#!/bin/bash
#infer_haplotype_freqs.sh 
### written by Susanne Tilk and Sharon Greenblum, Stanford University, 2018
#########################################
## NOTE: harp binary must be in path. download harp from https://bitbucket.org/dkessner/harp

usage()
{
    echo "usage: infer_haplotype_freqs.sh [ -b bamfile ] | [ -bl filename of bamfile list ]  
					  [ -s snps ] | [ -sl filename of snps list ]
					  [ -r refseq ]
					  [ -w wins=1000 ] 
					  [ -i illumina-quality-encoding?=0|1(default) ] 
					  [ -o outDir=<dirname bamfile> ] 
					  [ -t threads=1(default)] 
					  [ -h help ]
"
}
if [ $# -lt 1 ]; then usage; exit; fi

##### Main

wins=1000
illumina=1
threads=1

while [ "$1" != "" ]; do
    case $1 in
        -b | --bamfile )        shift
                                bams=($1)
                                ;;
        -bl | --bamlist )       shift
                                bams=($(cat $1))
                                ;;
        -s | --snps )           shift
                                snps=($1)
                                ;;
        -sl | --snpslist )      shift
                                snps=($(cat $1))
                                ;;
        -r | --refseq )         shift
                                refseq=$1
                                ;;
        -w | --wins )           shift
                                wins=$1
                                ;;
        -i | --illumina )       shift
                                illumina=$1
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

if [ -z "$bams" ] || [ -z "$snps" ] || [ -z "$refseq" ] ; then echo "Missing either bams or snps or refseq argument. Please try again"; usage; exit; fi
echo -e "********\ninferring haplotype freqs for [\n$(echo ${bams[*]} | tr ' ' '\n') \n]
using haplotypes in [\n${snps[*]}\n]
with ${wins}kb windows and $(echo $illumina | awk '($0==0){print "non-"}')illumina encoding\n*********"
if [ -z "$outDir" ]; then outDir="none"; fi

###################
## MAIN
#################################3

get_haps(){
	## VARS 
	wins=${1}
	illumina=${2}
	outDir=${3}
	refseq=${4}
	snp=${5}
	bam=${6}
	
	if [ "$outDir" == "none" ]; then outDir=$(dirname $bam); fi
	outFile=$outDir/$(basename $bam)
	
	## SET HARP WINDOWS
	likewindow=$(( $wins * 10000 ))  
	likestep=$(( $wins * 5000 )) 
	freqwindow=$(( $wins * 1000 )) 
	freqstep=$(( $wins * 100 )) 

	## GET CHROM START AND END POSITION 
	chrom=$(head -1 $snp | cut -f1 -d',')
	chrStart=1
	chrEnd=$(tail -n1 $snp | cut -d',' -f1)
	
	echo "running $bam vs. $snp"

	## RUN HARP IN EACH WINDOW
	for start in $(seq ${chrStart} ${likestep} ${chrEnd}) 
	do
	
	### get window definition
		(( stop=$start + $likewindow ))	
		if [ $stop -gt ${chrEnd} ]; then stop=${chrEnd}; fi
		
		### run harp like
		harp like \
		-b $bam \
		--refseq $refseq \
		--snps $snp \
		-r ${chrom}:${start}-${stop} \
		--stem ${outFile}.${chrom}_${start}_${stop} \
		$(echo $illumina | awk '($0==1){print "-I"}')  >/dev/null 2> /dev/null
		
		### run harp freq
		harp freq \
		-b $bam \
		--refseq $refseq \
		--snps $snp \
		-r ${chrom}:${start}-${stop} \
		--stem ${outFile}.${chrom}_${start}_${stop} \
		--window_step $freqstep \
		--window_width $freqwindow \
		$(echo $illumina | awk '($0==1){print "-I"}')  >/dev/null 2> /dev/null

		## CLEAN UP 
		rm -r ${outFile}.${chrom}_${start}_${stop}.output
		rm ${outFile}.${chrom}_${start}_${stop}.hlk
		
	done
	
	## CAT AND SORT HAPLOTYPE FREQUENCIES FROM ALL WINDOWS
	cat ${outFile}.${chrom}_*.freqs | tr ' ' '\t' | sort -k2g | tr '\t' ' ' > ${outFile}.${chrom}.freqs
	rm ${outFile}.${chrom}_*.freqs
	echo "harp freqs written to ${outFile}.${chrom}.freqs"
}

export -f get_haps
parallel --gnu -j${threads} get_haps $wins $illumina $outDir $refseq ::: ${snps[*]} ::: ${bams[*]}


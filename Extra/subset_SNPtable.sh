#!/bin/bash

###	Subsetting a .csv SNP file to selected columns/founderLines only
## GET PARAMS
usage()
{
    echo "usage: : subset_SNPtable.sh -s [snptable] -o [outfile] -f [founderfile] -k [keepallsites=false] -m [mincalls]"; 
    echo 'One or more variables are undefined. Enter the following: 
	1) snptable to be subsetted - first column is chrom position, second column is reference allele, rest of header row should list founder names  
	2) file name of resulting subsetted file 
	3) file name of list of founders (column headers) to keep in subsetted file 
	4) keep all sites regardless of whether they are segregating in the subset (default:false)
	5) minimum number of called genotypes to keep site'
}
if [ $# -lt 3 ]; then usage; exit; fi


##### Main

keep="false"

while [ "$1" != "" ]; do
    case $1 in
        -s | --snptable )       shift
                                snptable=($1)
                                ;;
        -o | --outfile )        shift
                                outfile=$1
                                ;;
        -f | --founderfile )    shift
                                founderfile=$1
                                ;;
        -k | --keepallsites )   shift
                                keep=$1
                                ;;
        -m | --mincalls )   	shift
                                mincalls=$1
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
## PRINT
echo "subsetting snp file $snptable to "$(cat $founderfile | wc -l )" selected founders and writing to $outfile"

## MAKE OUTDIR
outdir=$(dirname $outfile)
if [ ! -d ${outdir} ]; then
	mkdir -p ${outdir}
fi

## SUBSET
	## get field indices for selected lines ##
	fields="1,2,"$(head -1 $snptable | tr ',' '\n' | grep -nwf $founderfile | cut -f1 -d':' | tr '\n' ',' | sed 's/,$//')
	
	## extract and write to new file ##
	cat $snptable | cut -d',' -f$fields | awk -F',' -v keep="$keep" -v mincalls="$mincalls" '
		NR==1 {print $0}
		NR>1 { 
			if(keep!="false") {print $0}
			else {	
				ref=$2; 
				altCt=0; 
				refCt=0
				for (i=3;i<=NF;++i) {
					if (($i!=ref) && ($i!="N")) {altCt=altCt+1}
					if ($i==ref) {refCt=refCt+1}
				}
				if (altCt>0 && refCt>0 && (altCt+refCt)>=mincalls ) {print $0}
			}

		}' > $outfile

	### sanity check
	nfields_extract=$(( $(echo $fields | tr ',' ' ' | wc -w) - 2 ))
	nfields_new=$(( $(head -1 $outfile | awk -F "," '{print NF}') - 2 ))
	nsnps=$( tail -n +2 $outfile | wc -l )
	echo $nfields_extract " fields extracted, " $nfields_new " fields and $nsnps sites in new SNP file:
	$outfile"
	
	#####################



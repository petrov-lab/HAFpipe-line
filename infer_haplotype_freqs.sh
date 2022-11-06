#!/bin/bash
#infer_haplotype_freqs.sh
### written by Susanne Tilk and Sharon Greenblum, Stanford University, 2018
#########################################
## NOTE: harp binary must be in path. download harp from https://bitbucket.org/dkessner/harp

# ==================================================================================================
#      Usage
# ==================================================================================================

usage()
{
    echo "usage: infer_haplotype_freqs.sh
        [ -o --outdir ]     output directory for .${chrom}.freqs files;

        [ -d --scriptdir ]  directory in which HAF-pipe scripts are located;
                            #(default: directory of this script )

        [ -s --snptable ]   snp table to use for calculating haplotype and allele frequencies; tasks:2,3,4
                            #format:
                            #tab-delimited matrix of sites x founders, with 1 header line
                            #columns=position,ref-allele,founder1-allele,founder2-allele,...founderN-allele
                            #columnheaders=[nameofChrom](required!),Ref,[founder1-name],[founder2-name]...[founderN-name]
                            #nameofChrom must match a chromosome in the supplied reference fasta and a mapped chromosome in the supplied bamfile

        [ -m --method ]     method used for imputation;
                            #the script will look for the file ${snptable}.${method} and will use this to estimate haplotype frequencies
                            #if this file does not exist, an error will be thrown
                            #method can be '' (ie. to infer freqs with an un-imputed SNP table)

        [ -b --bamfile ]    name of bamfile with mapped reads for which haplotype frequencies will be inferred

        [ -r --refseq ]     reference fasta file

        [ -e --encoding ]   base quality encoding in bam file
                            #'illumina' (default)
                            #'sanger'

        [ -w --winsize ]    window size (in kb) for haplotype inference;
                            #(default: 1000)
                            #must be at least 10x smaller than length of chromosome

        [ -h help ]         show this help screen
"
}
if [ $# -lt 1 ]; then usage; exit; fi

# ==================================================================================================
#      Command Line Arguments
# ==================================================================================================

scriptdir=$(dirname "$0")
wins=1000
method=''
encoding="illumina"

while [ "$1" != "" ]; do
    case $1 in
        -d | --scriptdir )      shift
                                scriptdir=$1
                                ;;
        -o | --outdir )         shift
                                outdir=$1
                                ;;
        -s | --snptable )       shift
                                snptable=$1
                                ;;
        -m | --method )         shift
                                method="."$1
                                ;;
        -b | --bamfile )        shift
                                bamfile=$1
                                ;;
        -r | --refseq )         shift
                                refseq=$1
                                ;;
        -e | --encoding )       shift
                                encoding=$1
                                ;;
        -w | --winsize )        shift
                                wins=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     echo unknown flag $1 ; usage
                                exit 1
    esac
    shift
done

snptable=${snptable}${method}
if [ -z "$bamfile" ] || [ -z "$snptable" ] || [ -z "$refseq" ] ; then echo "Missing either bam or snptable or refseq argument. Please try again"; usage; exit; fi

echo -e "********\ninferring haplotype freqs for \n[ $bamfile ]\n
using haplotypes in \n[ $snptable ]\n
with ${wins}kb windows and $encoding base quality encoding\n*********"
if [ -z "$outdir" ]; then outdir=$(dirname $bamfile); fi
if [ ! -f ${snptable}.idx ]; then  $scriptdir/index_snp_table $snptable 50000; fi

# ==================================================================================================
#      Main
# ==================================================================================================

outfile=$outdir/$(basename $bamfile)

## SET HARP WINDOWS
likewindow=$(( $wins * 10000 ))
likestep=$(( $wins * 5000 ))
freqwindow=$(( $wins * 1000 ))
freqstep=$(( $wins * 100 ))

## GET CHROM START AND END POSITION
chrom=$(head -1 $snptable | cut -f1 -d',')
chrStart=1
chrEnd=$(tail -n1 $snptable | cut -d',' -f1)

## RUN HARP IN EACH WINDOW
for start in $(seq ${chrStart} ${likestep} ${chrEnd})
do

    ### get window definition
    (( stop=$start + $likewindow ))
    if [ $stop -gt ${chrEnd} ]; then stop=${chrEnd}; fi

    ### run harp like
    harp like \
    -b $bamfile \
    --refseq $refseq \
    --snps $snptable \
    -r ${chrom}:${start}-${stop} \
    --stem ${outfile}.${chrom}_${start}_${stop} \
    $(echo $encoding | awk '($0=="illumina"){print "-I"}')  >/dev/null 2> /dev/null

    ### run harp freq
    harp freq \
    -b $bamfile \
    --refseq $refseq \
    --snps $snptable \
    -r ${chrom}:${start}-${stop} \
    --stem ${outfile}.${chrom}_${start}_${stop} \
    --window_step $freqstep \
    --window_width $freqwindow \
    $(echo $encoding | awk '($0=="illumina"){print "-I"}')  >/dev/null 2> /dev/null

    ## CLEAN UP
    rm -r ${outfile}.${chrom}_${start}_${stop}.output
    rm ${outfile}.${chrom}_${start}_${stop}.hlk

done

## CAT AND SORT HAPLOTYPE FREQUENCIES FROM ALL WINDOWS
cat ${outfile}.${chrom}_*.freqs | tr ' ' '\t' | sort -k2g | tr '\t' ' ' > ${outfile}.${chrom}.freqs
rm ${outfile}.${chrom}_*.freqs
echo "harp freqs written to ${outfile}.${chrom}.freqs"

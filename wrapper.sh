### HAFpipe-line wrapper script
### written by Susanne Tilk and Sharon Greenblum, Stanford University, 2018
#########################################
## must be in path: 
#	1) index_snp_table
#	2) harp
#	3) parallel gnu (?)
#	4) tabix




usage()
{
    echo "usage: HAFpipe-line.sh [ -b bamfile ] | [ -bl filename of bamfile list ]  
				 [ -v vcf ]
				 [ -s steps to run ] 
				 [ -c chromosomes to run ]
				 [ -h help ]

    steps are:
	1 - make SNP table from VCF
	2 - impute SNP table
	3 - infer haplotype frequencies
	4 - calculate allele frequencies

     optional arguments:
	[ -f founder subset list: none(default)|filename ]	
	[ -e base quality encoding: illumina(default)|other ] 	
	[ -w window size (in kb) for haplotype inference: 10|100|1000(default) ]
	[ -o overwrite any existing files?: yes(default)|no ]
	[ -t nof threads/cores to use for parallization: 1(default) ]

     requirements:
	R (tested on versions >= 3.2)
	R libraries: 	
"
}


##### Main


while [ "$1" != "" ]; do
    case $1 in
        -b | --bamfile )        shift
                                bams=($1)
                                ;;
        -bl | --bamlist )       shift
                                bams=($(cat $1))
                                ;;
        -v | --vcf )            shift
                                vcf=$1
                                ;;
        -s | --steps )          shift
                                steps=($(echo $1 | tr ',' '\n' | sort )
                                ;;
        -c | --chroms )          shift
                                chroms=($(echo $1))
                                ;;
        -f | --founders )       shift
                                founders=$1
                                ;;
        -e | --encoding )       shift
                                enc=$1
                                ;;
        -w | --winsize )        shift
                                winsize=$1
                                ;;
        -o | --overwrite )      shift
                                overwrite=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done


# Confirm options 


###

for $step in ${steps[*]}; do
    case $step in
	1)	make_SNPtable_from_vcf.sh -v $vcf -c 





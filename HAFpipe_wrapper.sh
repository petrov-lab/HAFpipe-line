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
    echo "usage: HAFpipe-line.sh  [ -b bamfile ] OR [ -bl filename of bamfile list ] 
                    [ -r reference sequence ]
                    [ -v vcf ] OR [ -t snp table ]
                    [ -s steps to run (comma-separated)] 
                    [ -c chromosomes to run (comma-separated)]
                    [ -o output directory ]
                    [ -f founder subset list (default: none) ]
                    [ -m method to use for imputation (default:	'simpute'	ie. simple impututation
							alt:	'npute' 	see Roberts et al., 2007 - doi:10.1093/bioinformatics/btm220 
		    ]
		    [ -n number of neighboring sites to use with 'npute' imputation (default: 20)]
		    [ -e base quality encoding (default: illumina) ]
                    [ -w window size (in kb) for haplotype inference (default: 1000kb) ]
                    [ -p number of threads to use for parallization  (default: 1) ]
                    [ -h help ]
    steps are:
	1 - make SNP table from VCF
	2 - subset SNP table
	3 - impute SNP table
	4 - infer haplotype frequencies
	5 - calculate allele frequencies


     requirements:
	R (tested on versions >= 3.2)
	R libraries: 	
"
}


##### Main


## Set Default Parameters
founders="all"
method="simple"
nSites=20
enc="-I"
winsize=1000
threads=1
subsetName=""


## Parse User Parameters
while [ "$1" != "" ]; do
    case $1 in
        -b | --bamfile )        shift
                                bamstring="-b $1"
                                ;;
        -bl | --bamlist )       shift
                                bamstring="-bl $1"
                                ;;
        -r | --refseq )         shift
                                ref=$1
                                ;;
        -v | --vcf )            shift
                                vcf=$1
                                ;;
        -t | --snptable )       shift
                                snptable=$1
                                ;;
        -s | --steps )          shift
                                steps=($(echo $1 | tr ',' '\n' | sort ))
                                ;;
        -c | --chroms )          shift
                                chroms=($(echo $1))
                                ;;
     	-o | --outdir )         shift
                                outDir=$1
                                ;;
        -fl | --founderList )   shift
                                flist=$1
                                ;;
        -f | --founderSetName )	shift
                                fset="."$1
                                ;;
        -m | --method )         shift
                                method=$1
                                ;;
        -n | --nSites )         shift
                                nSites=$1
                                ;;
        -e | --encoding )       shift
                                enc=$1
                                ;;
        -w | --winsize )        shift
                                winsize=$1
                                ;;
        -p | --par )            shift
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


# Confirm options 
## TO DO!


### BEGIN
for chrom in ${chroms[*]}; do
	for step in ${steps[*]}; do
    		case $step in
		1)	snptable=$(make_SNPtable_from_vcf.sh -v $vcf -c $chrom -o $outDir | tail -1) 
			#optional args: [-f firstSampleCol=10] [-m minNofCalls=2]
			;;  
		2)	subset_SNPtable.sh ${snptable} ${snptable}${fset} $flist
			;;
		3)	case $method in
			simpute) impute_SNPtable.sh ${snptable}${fset}; 
				;;
			npute)	npute_SNPtable.sh ${snptable}${fset} $nSites
				;;
			* )	echo "not a valid imputation method"
				exit 1
			esac
			;;
		4)	infer_haplotype_freqs.sh $bamstring -s ${snptable}${fset}.${method} -r $refseq -w $wins -i $enc -o $outDir -t $threads
			;;
		5)	calculate_HAFs.sh $bamstring -s ${snptable}${fset} -o $outDir -t $threads
                        ;;
        	* )      usage
                        exit 1
    		esac
	done
done





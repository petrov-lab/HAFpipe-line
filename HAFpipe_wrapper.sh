### HAFpipe-line wrapper script
### written by Susanne Tilk and Sharon Greenblum, Stanford University, 2018
#########################################


usage()
{
    echo "usage: HAFpipe-line.sh  \
        [ -t --tasks ]      tasks to run (comma-separated)

        [ -l --logfile ]    name of file to write log of commands to

        [ -d --scriptdir ]  directory in which HAF-pipe scripts are located; tasks:1,2,3,4
                            #(default: directory of this script )

        [ -o --outdir ]     output directory; tasks:1,3,4 

        [ -v --vcf ]        vcf file to be converted to snp table; tasks:1

        [ -c --chrom ]      name of chromosome to extract from vcf; tasks:1

        [ -k --keephets ]   whether to keep heterozygous calls as ambiguous bases in SNPtable (rather than treat them as missing and impute them)

        [ -s --snptable ]   snp table to use for calculating haplotype and allele frequencies; tasks:2,3,4 
                            #will be overwritten if task 1 is run in conjuction with other tasks

        [ -m --method ]     method to use for imputation in task 2 or file extension for task 3; tasks:2,3
        					for task 2, method must be one of:
                            #'simpute' (simple imputation)
                            #'npute' (see Roberts et al., 2007 - doi:10.1093/bioinformatics/btm220 )
                            for task 3, method can be any string:
                            #default:'none' 
                            #if a string other than 'none' is supplied, the script will look for the file called [snptable].[method] and will use this to infer haplotype frequencies
                            #if method is not specified or is 'none', the original [snptable] with potential missing calls will be used to infer haplotype frequencies

        [ -n --nsites ]     number of neighboring sites to use with 'npute' imputation; tasks:2
                            #(default: 20)

        [ -b --bamfile ]    name of bamfile with mapped reads; tasks:3,4 

        [ -r --refseq ]     reference sequence; tasks:3

        [ -e --encoding ]   base quality encoding in bam files; tasks:3
                            #'illumina' (default)
                            #'sanger'

	[ -g --generations ] number of generations of recombination; used to calculate window size for haplotype inference; tasks=3
	
	[ -a --recombrate ] recombination rate used to calculate window size for haplotype inference; tasks=3	

	[ -q --quantile ]   quantile of expected unrecombined segment distribution to use for determining haplotype inference window size; tasks:3

        
	[ -w --winsize ]    user-defined window size (in kb) for haplotype inference; tasks:3
                            #(overrides -g and -a) 

        [ -h help ]         show this help screen
    
    #tasks are:
	1 - make SNP table from VCF
	2 - impute SNP table
	3 - infer haplotype frequencies
	4 - calculate allele frequencies

    #requirements:
	R (tested on versions >= 3.2)
	R libraries: (data.table)
    for imputation method 'npute':
    -- python 2
    -- numpy
    
    ## must be installed and on path: 
   1) harp
   2) tabix
   3) bgzip
"
}


##### Main


## Set Default Parameters
tasks=0
logfile=$(echo HAFpipe-log.`date +%Y-%m-%d.%H%M%S`)
scriptdir=$(dirname "$0")
outdir=.
vcf=""
chrom=""
keephets=""
method=""
snptable=""
nsites=20
bamfile=""
refseq=""
encoding="illumina"
winsize=""
recombrate=.0000000239
quantile=18
gens=""



if [ "$1" == "" ]; then usage; exit; fi

## Parse User Parameters
while [ "$1" != "" ]; do
    case $1 in
        -t | --tasks )          shift
                                tasks=($(echo $1 | tr ',' '\n' | sort ))
                                ;;
        -l | --logfile )        shift
                                logfile=$1
                                ;;
        -d | --scriptdir )      shift
                                scriptdir=$1
                                ;;
        -o | --outdir )         shift
                                outdir=$1
                                ;;
        -v | --vcf )            shift
                                vcf=$1
                                ;;
        -c | --chrom )          shift
                                chrom=$1
                                ;;
        -k | --keephets )       
                                keephets="--keephets"
                                ;;
        -s | --snptable )       shift
                                snptable=$1
                                ;;
        -m | --method )         shift
                                case $1 in 
                                    none ) method="" ;;
				                    * ) method="."$1 ;;
				                esac
                                ;;
        -n | --nsites )         shift
                                nsites=$1
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
        -g | --generations )    shift
                                gens=$1
                                ;;
        -a | --recombrate )     shift
                                recombrate=$1
                                ;;
        -q | --quantile )	    shift
                                quantile=$1
                                ;;
        -w | --winsize )        shift
                                winsize=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     echo unknown flag $1 ; usage
                                exit 1
    esac
    shift
done


# Write options to log file
echo "
        $(date)
        ####PARAMETERS########
        --tasks: $tasks
        --scriptdir $scriptdir
        --outdir $outdir
        --vcf $vcf
        --chrom $chrom
        --snptable $snptable
        --method $method
        --nsites $nsites
        --bamfile $bamfile
        --refseq $refseq
        --encoding $encoding
        --generations $gens
        --recombrate $recombrate
        --quantile $quantile
        --winsize $winsize

        #####COMMANDS#######



" > $logfile; 


### BEGIN
for task in ${tasks[*]}; do
	case $task in
	1)	echo -e "$scriptdir/make_SNPtable_from_vcf.sh -v $vcf -c $chrom -o $outdir $keephets \n $scriptdir/count_SNPtable.sh $snptable \n $scriptdir/numeric_SNPtable.sh $snptable" >> $logfile
		snptable=$($scriptdir/make_SNPtable_from_vcf.sh -v $vcf -c $chrom -o $outdir $keephets| tail -1) #optional args: [-f firstSampleCol=10] [-m minNofCalls=2]
		$scriptdir/count_SNPtable.sh $snptable
		Rscript $scriptdir/numeric_SNPtable.R $snptable
		;;  
	2)	case $method in
		".simpute") echo "$scriptdir/impute_SNPtable.sh ${snptable}" >> $logfile; $scriptdir/impute_SNPtable.sh ${snptable} 
		;;
		".npute") echo "$scriptdir/npute_SNPtable.sh ${snptable} $nsites" >> $logfile; $scriptdir/npute_SNPtable.sh ${snptable} $nsites
		;;
		* )	echo "not a valid imputation method"
		exit 1
		esac
		;;
    	3)	if [ "$winsize" == "" ]; then
			if [ "$gens" != "" ]; then
				if [ $quantile -lt 100 ] && [ $quantile -gt 0 ]; then
					chromlength=$(tail -1 $snptable | cut -f1 -d',')
					winsize=$(Rscript -e \
						"args <- commandArgs(TRUE);args=as.numeric(args); R=args[1];L=args[2];G=args[3];Q=args[4]; round(qexp(Q/100,1/(L/(R*L*G+1)))/1000)" \
				$recombrate $chromlength $gens $quantile | cut -f2 -d' ');
                			echo "running harp with window size: $winsize kb"
				else echo "invalid quantile value. must be between 0 and 100."; exit; 
				fi
			else echo "number of generations must be defined to calculate window size"; exit
			fi	
		fi
		echo "$scriptdir/infer_haplotype_freqs.sh -b $bamfile -s ${snptable}${method} -r $refseq -w $winsize -e $encoding -o $outdir -d $scriptdir" >> $logfile
		$scriptdir/infer_haplotype_freqs.sh -b $bamfile -s ${snptable}${method} -r $refseq -w $winsize -e $encoding -o $outdir -d $scriptdir >> $logfile
		;;
	4)	chrom=$(head -1 $snptable | cut -f1 -d',')
		freqs=$outdir/$(basename $bamfile)".$chrom.freqs"
		if [ ! -e $freqs ]; then echo "$freqs not found!"; exit 1; fi
		echo "$scriptdir/calculate_HAFs.sh -f $freqs -s $snptable -o $outdir" >> $logfile
        $scriptdir/calculate_HAFs.sh -f $freqs -s $snptable -o $outdir -d $scriptdir >> $logfile
        ;;
    * ) usage
        exit 1
    esac
done





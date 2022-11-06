    #!/bin/bash
#####
# SIMULATE POOLED SEQUENCING READS FROM A FORQS RUN - one chrom only
########


### SET DEFAULTS
outDir=.
poolSize=100
errRate=.002
maindir=$(dirname "$0")/..

# ==================================================================================================
#      Usage
# ==================================================================================================

usage()
{
    echo "
    usage: HAFpipe-sim.run_forqs.sh
    -f | --forqsDir)
    -p | --pop )
    -n | --poolSize )
    -g | --gen )
    -c | --covList )
    -s | --snps )
    -f | --refFasta )
    -e | --errRate )
    -o | --out )
    -h | --help )
"
}

# ==================================================================================================
#      Command Line Arguments
# ==================================================================================================

if [ "$1" == "" ]; then
    usage
    exit 1
fi

## Parse User Parameters
while [ "$1" != "" ]; do
    case $1 in
        -f | --forqsDir )       shift
                                forqsDir=$1
                                ;;
        -p | --pop )            shift
                                pop=$1
                                ;;
        -n | --poolSize )       shift
                                poolSize=$1
                                ;;
        -g | --gen )            shift
                                gen=$1
                                ;;
        -c | --covList )        covString=$1
                                ;;
        -s | --snps )           shift
                                snps=$1
                                ;;
        -r | --refFasta )       shift
                                refFasta=$1
                                ;;
        -e | --errRate )        shift
                                errRate=$1
                                ;;
       -o | --out )             shift
                                outDir=$1
                                ;;
       -h | --help )           usage
                                exit 1
                                ;;
        * )                     echo unknown flag $1 ; usage
                                exit 1
    esac
    shift
done

# ==================================================================================================
#      Main
# ==================================================================================================

chrom=$(head -1 $snps | cut -f1 -d',')
runID=$(basename $forqsDir | sed 's/forqs.//' | cut -f1 -d'_')
round=$(basename $forqsDir | sed 's/forqs.//' | cut -f2 -d'_')
snpsDir=$outDir/$runID/snpTables
bamDir=$outDir/$runID/bams
logfile=$outDir/logs/run_${runID}_sim.log

####################
## recombine snp table to make sampled chroms and simulate reads from it
#####################

## recombine snp table
if [ ! -e  $snpsDir/round${round} ] ; then mkdir -p $snpsDir/round${round}; fi

echo "recombining snps in
$snps
for $runID gen ${gen} pop ${rep} round ${round}" >> $logfile #> /dev/tty #

${maindir}/simulations/recombineSNPtable_reps.sh $snps $forqsDir $gen $rep $poolSize $snpsDir/round${round}
recombinedSNPs=$snpsDir/round${round}/$(basename $snps).g${gen}.r${rep}
if [ -e ${recombinedSNPs}.idx ]; then rm ${recombinedSNPs}.idx; fi

${maindir}/scripts/index_snp_table $recombinedSNPs 10000
${maindir}/scripts/count_SNPtable.sh $recombinedSNPs
mv $recombinedSNPs.alleleCts $snpsDir/$(basename $snps).g${gen}.r${rep}.round${round}.alleleCts
mv $recombinedSNPs.brkpts $snpsDir/$(basename $snps).g${gen}.r${rep}.round${round}.brkpts
echo "true allele frequencies calculated: " >> $logfile #> /dev/tty
ls -lh $snpsDir/$(basename $snps).g${gen}.r${rep}.round${round}.alleleCts >> $logfile

## simulate reads
if [ ! -e  $bamDir ] ; then mkdir -p $bamDir; fi
for cov in $(echo $covString | tr ',' '\n'); do
    simreadsStem=$bamDir/simreads.${gen}g.r${rep}.${cov}x.round${round}
    echo "simulating reads at ${cov}x cov" >> $logfile
    ${maindir}/simulations/runSimreads.sh $recombinedSNPs $simreadsStem $cov $errRate $refFasta
    echo "reads written to:" >> $logfile
    ls -lh $simreadsStem.bam >> $logfile
done
rm $recombinedSNPs

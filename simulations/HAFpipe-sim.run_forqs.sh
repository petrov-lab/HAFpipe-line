#!/bin/bash
#####
# RUN FORQS RECOMBINATION SIMULATION - one chrom only
########

### SET DEFAULTS
runID=default
rounds=1
pops=1
popSize=1000
recMap=dmel_recRates_2L-100kb.csv
gens=200
loci=0
selection=0
constantEffect=1
report_every_nGens=1
snps=99.clean.SNP.HARP.segregating.imputedSim.numeric
outDir="."

# ==================================================================================================
#      Usage
# ==================================================================================================

usage()
{
    echo "
    usage: HAFpipe-sim.run_forqs.sh
     -i | --runID )
    -r | --rounds )
    -p | --pops )
    -n | --popSize )
    -m | --recMap )
    -g | --gens )
    -l | --loci )
    -c | --selectionCoeff )
    -s | --snps )
    -v | --report )
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
        -i | --runID )          shift
                                runID=$1
                                ;;
        -r | --rounds )         shift
                                rounds=$1
                                ;;
        -p | --pops )           shift
                                pops=$1
                                ;;
        -n | --popSize )         shift
                                popSize=$1
                                ;;
        -m | --recMap )         shift
                                recMap=$1
                                ;;
        -g | --gens )          shift
                                gens=$1
                                ;;
        -l | --loci )           shift
                                loci=$1
                                ;;
        -c | --selectionCoeff ) selection=$1
                                ;;
        -s | --snps )           shift
                                snps=$1
                                ;;
        -v | --report )         shift
                                report_every_nGens=$1
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

### RUNS FROM PAPER
outDir=forqs/

#runID=sus99_neutral; pops=5; loci=0; selection=0; gens=50; rounds=3
runName=sus99_long; pops=3; loci=5; selection=.025; gens=250; rounds=1; report_every_nGens=25
#runID=sus99_10sites_strong; pops=3; loci=10; selection=.1; gens=50; rounds=5
#runID=sus99_5sites_strong; pops=5; loci=5; selection=.1; gens=50; rounds=3
#runID=sus99_5sites_weak; pops=5; loci=5; selection=.025; gens=50; rounds=3
#runID=sus99_10sites_weak; pops=5; loci=10; selection=.025; gens=50; rounds=3
#runID=dgrp_5sites_weak; pops=5; loci=5; selection=.025; gens=50; rounds=3; snps=freeze2.2L.snptable.imputed.numeric
#runID=cendr_5sites_weak; pops=5; loci=5; selection=.025; gens=50; rounds=3; snps=CeNDR_100.chrI.snptable.imputed.numeric; recMap=CeNDR_recombinationMap_chrI.txt

# ==================================================================================================
#      Main
# ==================================================================================================

if [ ! -e $outDir ]; then mkdir $outDir; fi
for round in `seq 1 $rounds`; do
seed=$round

configFile=$outDir/forqs.${runID}_${seed}.config.txt

echo "
creating forqs config file for:
$pops populations of size $popSize,
founded with $(head -1 $snps | tr ',' '\t' | cut -f2- | wc -w) haplotypes from $snps,
each with $loci selected loci with strength $selection,
recombining for $gens generations,
according to rec rates in $recMap.
writing to $outDir"

###############
##  SET UP
##############
chromLength=$(tail -1 $recMap | tr ' ' '\t' | cut -f1)

echo "
# forqs config file

# set up population
PopulationConfigGenerator_ConstantSize pcg
    population_size = $popSize
    generation_count = $gens
    population_count = $pops
    chromosome_pair_count = 1
    chromosome_lengths = $chromLength " > $configFile
    if [ $loci -gt 0 ]; then
    echo "    fitness_function = fitness " >> $configFile ##
    fi
echo "
###############
# set recomb map
###############
# generate recombination positions based on a genetic map
RecombinationPositionGenerator_RecombinationMap rpg_map
    filename = $recMap

  # report population haplotypes after each set of n generations
    Reporter_Population reporter_population
    update_step = $report_every_nGens
" >> $configFile

##############
## ADD SELECTED LOCI
###############
if [ $loci -gt 0 ]; then

## create ms file to:
## assign loci to positions in genome
## assign alleles to haplotypes
msfile=$outDir/forqs.${runID}_${seed}.ms.txt
tail -n +2 $snps | shuf | head -$loci | tr ',' '\t' | sort -k1g > alleles.tmp
echo "
//
segsites: $loci
positions: "$(cat alleles.tmp | cut -f1 -d',' |  awk -v chromLength="$chromLength" '{printf "%f ",$ii/chromLength}') > $msfile

cat alleles.tmp | awk -v popSize="$popSize" -v nPops="$pops" '
{       # collect alleles for each founder at selected sites
    for(ii=2;ii<=NF;ii++){ hapAlleles[ii-1]=hapAlleles[ii-1]$ii }
}
END { # print allele string for each forqs index
    for(forqsIX=0;forqsIX<(2*popSize*nPops);forqsIX++){
        founderIX=1+((forqsIX-forqsIX%2)/2)%(NF-1)
        print hapAlleles[founderIX]
    }
}
' >> $msfile

## add loci to forqs config file
echo "
LocusList selected_loci " >> $configFile
cat alleles.tmp | cut -f1 -d',' | awk '{print "chromosome:position = 1 "$1}' >> $configFile

echo "
VariantIndicator_File variants
    msfile = $msfile
    loci = selected_loci
" >> $configFile

## set selection coefficients
if [ $constantEffect -gt 0 ]; then echo "
Distribution_Constant positive_effect
    value = $selection
" >> $configFile ; else echo "
do something else here for non-uniform selection coeffs across sites?
" >> $configFile
fi

## set up quantitative trait with independent loci effects
echo "

Distribution_Constant no_dominance
    value = 0

QTLEffectGenerator generator_positive
    locus_list = selected_loci
    effect_size_distribution = positive_effect
    dominance_distribution = no_dominance

QuantitativeTrait_IndependentLoci qt
    qtl_effect_generator = generator_positive
    environmental_variance = 0.01


# determine how next gen is selected (ie selection)
# OPT USED: optimum - assign fitness 1 to optimum trait value, with fitness decaying linearly away from the optimum
# other opts: truncation - assign fitness 1 if trait value exceeds threshold, 0 otherwise. Threshold is determined each generation by the user-specified proportion of individuals to be selected to reproduce.
FitnessFunction_Optimum fitness
    quantitative_trait = qt
    optimum = 1
    radius:power = 1 1

#report output from simulation
Reporter_AlleleFrequencies reporter_allele_freqs
    locus_list = selected_loci

Reporter_TraitValues reporter_trait_values
    quantitative_traits = qt
"  >> $configFile
fi

############
## PUT IT ALL TOGETHER
###########
# SimulatorConfig is the top-level module, and must appear last in the
# configuration file.  This module specifies which of the previously
# specified primary modules are plugged into the main simulator.

echo "
SimulatorConfig
    output_directory = $outDir/forqs.${runID}_${seed}
    population_config_generator = pcg
    write_popconfig = 1 # set this flag to write out the popconfig file
    recombination_position_generator = rpg_map
    reporter = reporter_population
    seed = $seed"  >> $configFile
if [ $loci -gt 0 ]; then echo "    quantitative_trait = qt
    quantitative_trait = fitness
    reporter = reporter_allele_freqs
    reporter = reporter_trait_values
    variant_indicator = variants
" >> $configFile
fi

forqs $configFile > forqs.log

echo "-->results for round $round written to
$outDir/forqs.${runName}_${seed}"
done

#!/bin/bash


## VARS
#########
snps=${1}
outFile=${2}
cov=${3}
err=${4}
ref=${5}
indivProportionsFile=${6:-NULL}. #ie. /mnt/lab_NGS_backup_v2/poolTest2/04_harp_Runs/HARP.SNP.input/final.universe/singletons.all.chrom.SNP.HARP.founder_identifiers.sharonProportions
nondros=${7:-0}
simreads=${8:-/mnt/cages/sim/scripts/simreads/simreads}

outDir=$(dirname $outFile)
outBase=$(basename $outFile)

########
# SET UP PROPORTION OF CHROMOSOME EACH IN POOL
#######
chrom=$(head -1 $snps | tr ',' '\t' | cut -f1 )
nofFounderChroms=$(( $(head -1 $snps | tr ',' ' ' | wc -w) - 2 ))
if [ -e $indivProportionsFile ]; then
        echo "pooling individuals in SNP table according to proportions in $indivProportionsFile"
        haplotypes=$(echo $(cat $indivProportionsFile | cut -f2 -d' '))
    else
        echo "pooling $nofFounderChroms evenly for sequencing"
        haplotypes=$(yes "$(echo $nofFounderChroms | awk '{print 1/$1}') " | head -n $nofFounderChroms | tr -d '\n') ## EVEN POOL
fi
ref_end=$(samtools faidx $ref $chrom | tail -n +2 | sed 's/\S//' |  wc -c)
regionEnd=$(( $ref_end - 451 ))  ## can only simulate with enough room to make a paired end plus insert

#bqiFile=/mnt/cages/scripts/simreads_out/R1.8F.bqi
#filename_bqi $bqiFile
#region $chrom:1-$regionEnd

echo "

filename_refseq ${ref}

filename_snps ${snps}

filename_stem $outBase

region $chrom:150-$regionEnd

haplotype_frequencies $haplotypes

coverage $cov

error_rate $err

read_length 150

" > $outFile.config.txt
################
# RUN SIMREADS [RAW]
##########
dir=$(pwd)
cd $outDir
$simreads $outFile.config.txt
cd $dir


###########
# MOVE TO OUTDIR AND INDEX
##########

if [ $nondros > 0 ]; then
    echo "@SQ	SN:$chrom	LN:$ref_end" > ${outFile}.sam.new
    echo "@PG	ID:simreads PN:simreads VN:1.0" >> ${outFile}.sam.new
    grep -v "^@" ${outFile}.sam >> ${outFile}.sam.new

    mv ${outFile}.sam.new ${outFile}.sam
fi

samtools view -Suh ${outFile}.sam | \
samtools sort -@10 -o ${outFile}.bam -
samtools index ${outFile}.bam
rm ${outFile}.sam
rm ${outFile}.actual.freqs
rm ${outFile}.true.freqs
rm ${outFile}.seed

#!/bin/bash


## VARS
#########
snps=${1}
outFile=${2} ## MUST SIMULATE TO HOME DIR, AND THEM MOVE OUTPUT TO SPECIFIED DIR
outDir=$(dirname $outFile)
outBase=$(basename $outFile)
cov=${3}
err=${4}
ref=/mnt/Alan/dmel_reference/all_dmel.fasta

########
# SET UP PROPORTION OF CHROMOSOME EACH IN POOL
#######
chrom=$(head -1 $snps | tr ',' '\t' | cut -f1 )
nofFounderChroms=$(( $(head -1 $snps | tr ',' ' ' | wc -w) - 2 ))
haplotypes=$(yes "$(echo $nofFounderChroms | awk '{print 1/$1}') " | head -n $nofFounderChroms | tr -d '\n') ## EVEN POOL
dmelRef_end=$(cat $ref | grep ">$chrom" | awk '{split($3,loc,".");print loc[3]}' | sed 's/;//')
regionEnd=$(( $dmelRef_end - 451 ))  ## can only simulate with enough room to make a paired end plus insert

#bqiFile=/mnt/cages/scripts/simreads_out/R1.8F.bqi
#filename_bqi $bqiFile

echo "

filename_refseq ${ref}

filename_snps ${snps}

filename_stem $outBase

region $chrom:1-$regionEnd
haplotype_frequencies $haplotypes

coverage $cov

error_rate $err

read_length 150

" > $outFile.config.txt
################
# RUN SIMREADS [RAW]
##########
/mnt/cages/scripts/simreads $outFile.config.txt


###########
# MOVE TO OUTDIR AND INDEX
##########
mv ${outBase}.sam $outDir
samtools view -Suh ${outFile}.sam | \
samtools sort -@10 -o ${outFile}.bam -
rm ${outFile}.sam
samtools index ${outFile}.bam
rm ${outBase}*


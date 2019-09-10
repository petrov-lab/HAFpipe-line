#!/bin/bash

## SET PARAMETERS
threads=10
doSUS99=1
doDGRP=0
doCENDR=0
doMakeBams=1
doCalcHAFs=1

rootDir=/mnt/cages/sus99/selection
if [ ! -e  $rootDir/logs ] ; then mkdir -p $rootDir/logs; fi
logfile=$rootDir/logs/sim_harp_af.selection.log
runID=$RANDOM
echo "#### run$runID started "`date +%Y-%m-%d` >> $logfile

#### SUS99
if [ "$doSUS99" -gt 0 ]; then
	echo "running SUS99"
	#regimes=( sus99_5sites_weak sus99_5sites_strong sus99_10sites_weak sus99_long);
	regimes=( sus99_long );
	reps=(1 2 3); # reps=(1 2 3 4 5)
	rounds=(1) #(1 2 3)
	gens=(25 50 75 100 125 150 175 200 225 250) #( 5 10 15 25 50);
	covs=(5) #(1 3 5 7 9)
	pctMissings=(0) # pctMissings=(0 1 2 3 5 10)
	snps=/mnt/lab_NGS_backup_v2/poolTest2/04_harp_Runs/HARP.SNP.input/final.universe/99.clean.SNP.HARP.segregating.imputedSelectionSim
	refFasta=/mnt/cages/ref_data/d_mel/fasta_all/all_dmel.fasta
	recRate=.0000000239

	# so far: 
	## bams: 450 rep1 round1-3 genall cov1-10 = 150  (extra: 75 1g files 2/3) created 3/31/2019
	## 		     rep2-5 round1-3 genAll cov1,3,5,7,9 = 300 created 2/3/2019
	## 0 missing: total 210
	##		rep1 round1-3 genAll cov1-10 = 150 created 3/31
	## 		rep2-5 round1-3 genAll cov5 = 60 created 2/7
	## 1 missing: total 175 
	## 		rep1 round1-3 genAll cov1,3,5,7,9 = 75 created 3/31/2019
	##      rep2-5 round1 genAll cov1,3,5,7,9 (extra: 20 10x files 2/7) = 100 created 2/10/2019 and 5/8-9/2019 
	## 2 missing: total 175
	## 		rep1 round1-3 genAll cov1,3,5,7,9 = 75 created 3/31
	##      rep2-5 round1 genAll cov1,3,5,7,9 = 100 created 5/8-9
	## 3 missing: total 175
	## 		rep1 round1-3 genAll cov1,3,5,7,9 = 75 created 3/31
	##      rep2-5 round1 genAll cov1,3,5,7,9 = 100 created 5/8-9
	## 5 missing: total 175
	## 		rep1 round1-3 genAll cov1,3,5,7,9 = 75 created 3/31-4/1
	##      rep2-5 round1 genAll cov1,3,5,7,9 (extra: 20 10x files 2/7) = 70 created 2/10/2019 and 5/8-9/2019 
	## RUNNING rep2-5 covs3,7,9
	## 10 missing: total 115
	## 		rep1 round1-3 genAll cov1,3,5,7,9 = 75 created 4/1
	##      rep2-5 round1 genAll cov1,5(,10) = 40
	## RUNNING rep2-5 covs3,7,9
	
	##RUN rep2-5 round1 genAll cov3,7,9 pct1,2,3,4,10 = 300
fi


#### DGRP
if [ $doDGRP -gt 0 ]; then
	echo "running DGRP"
	regimes=(dgrp_5sites_weak) 
	reps=(1)
	rounds=(2); #(1 2)
	gens=(5) #(5 15 50); 
	covs=(1 3 5) #(1 3 5)
	pctMissings=(1 5 10)
	snps=/mnt/cages/ref_data/dgrp/snps/all_lines/freeze2.2L.snptable.imputed
	refFasta=/mnt/cages/ref_data/d_mel/fasta_all/all_dmel.fasta
	recRate=.0000000239
fi

### CENDR
if [ $doCENDR -gt 0 ]; then
	echo "running CENDR"
	regimes=(cendr_5sites_weak) 
	reps=(1)
	rounds=(1 2 3)
	gens=(5 15 50)  #gens=(5 10 15 25 50); 
	covs=(1 3 5)
	pctMissings=(1 5 10) #pctMissings=(0 1 2)
	snps=/mnt/lab_NGS_backup_v2/poolTest2/10_CeNDR/snpTables/CeNDR.chrI.clean.snptable_100_totalSubsetted.segregating.forMakingReads.imputed
	refFasta=/mnt/lab_NGS_backup_v2/poolTest2/10_CeNDR/fastas/Caenorhabditis_elegans.WBcel235.dna.chromosome.I.fa
	recRate=.0000000333 
fi

covString=$(echo $(for cov in ${covs[*]}; do echo $cov; done) | tr ' ' ',');

makeBams() {
	rootDir=/mnt/cages/sus99/selection
	regime=$1
	round=$2
	rep=$3
	gen=$4
	covString=$5
	snps=$6
	refFasta=$7
	logfile=$rootDir/logs/run_${regime}_sim.log
	
	poolSize=100
	errRate=.002
	scriptDir=/mnt/cages/scripts/HAFpipe-line/
	chrom=$(head -1 $snps | cut -f1 -d',')

	snpsDir=$rootDir/$regime/snpTables
	bamDir=$rootDir/$regime/bams


	####################
	## recombine snp table to make sampled chroms and simulate reads from it
	#####################
	## RUN FORQS BEFORE RUNNING THIS SCRIPT!
	forqsDir=$rootDir/forqs/forqs.${regime}_${round}
	
	## recombine snp table 
	if [ ! -e  $snpsDir/round${round} ] ; then mkdir -p $snpsDir/round${round}; fi
	echo "recombining snps in 
	$snps
	for $regime gen ${gen} rep ${rep} round ${round}" >> $logfile #> /dev/tty #
	$scriptDir/simulations/recombineSNPtable_reps.sh $snps $forqsDir $gen $rep $poolSize $snpsDir/round${round}
	recombinedSNPs=$snpsDir/round${round}/$(basename $snps).g${gen}.r${rep}
	if [ -e ${recombinedSNPs}.idx ]; then rm ${recombinedSNPs}.idx; fi
	$scriptDir/index_snp_table $recombinedSNPs 10000
	$scriptDir/count_SNPtable.sh $recombinedSNPs
	mv $recombinedSNPs.alleleCts $snpsDir/$(basename $snps).g${gen}.r${rep}.round${round}.alleleCts
	mv $recombinedSNPs.brkpts $snpsDir/$(basename $snps).g${gen}.r${rep}.round${round}.brkpts
	echo "true allele frequencies calculated: " >> $logfile #> /dev/tty
	ls -lh $snpsDir/$(basename $snps).g${gen}.r${rep}.round${round}.alleleCts > /dev/tty #>> $logfile
	
	## simulate reads
	if [ ! -e  $bamDir ] ; then mkdir -p $bamDir; fi
	for cov in $(echo $covString | tr ',' '\n'); do
		simreadsStem=$bamDir/simreads.${gen}g.r${rep}.${cov}x.round${round}
		echo "simulating reads at ${cov}x cov" > /dev/tty #>> $logfile
		$scriptDir/simulations/runSimreads.sh $recombinedSNPs $simreadsStem $cov $errRate $refFasta
		echo "reads written to:" > /dev/tty #>> $logfile
		ls -lh $simreadsStem.bam > /dev/tty #>> $logfile
	done	
	rm $recombinedSNPs
}
export -f makeBams
if [ $doMakeBams -gt 0 ]; then
	echo "run$runID BAMS: parallel --gnu -j$threads makeBams ::: ${regimes[*]} ::: ${rounds[*]} ::: ${reps[*]} ::: ${gens[*]} ::: $covString ::: $snps ::: $refFasta" >> $logfile
	parallel --gnu -j$threads makeBams ::: ${regimes[*]} ::: ${rounds[*]} ::: ${reps[*]} ::: ${gens[*]} ::: $covString ::: $snps ::: $refFasta
	echo "run$runID BAMS: complete" >> $logfile
fi

calcHAFs() {
	#########################
	## VARS
	#################
	regime=$1
	round=$2
	rep=$3
	gen=$4
	cov=$5
	pctMissing=$6;
	snps=$7
	refFasta=$8
	recRate=$9
	poolSize=100
	encoding=illumina
	chrom=$(head -1 $snps | cut -f1 -d',')
	quantiles=(18)
	

	#########################
	## DIRS	
	#################
	scriptDir=/mnt/cages/scripts/HAFpipe-line/
	rootDir=/mnt/cages/sus99/selection
	logfile=$rootDir/logs/run_${regime}_sim.log #
	snpsDir=$rootDir/$regime/snpTables
	bamDir=$rootDir/$regime/bams
	
	for quantile in ${quantiles[*]}; do
		
	HAFDir=$rootDir/$regime/HAFs_${pctMissing}missing_q${quantile}
	if [ ! -e  $HAFDir ] ; then mkdir -p $HAFDir; fi
	
	#########################
	## run harp and calculate HAFs with and without missing sites
	#################
	simreadsStem=$bamDir/simreads.${gen}g.r${rep}.${cov}x.round${round}
	
	echo "calculating HAFs for bam file:" >> $logfile 
	ls -lh ${simreadsStem}.bam >> $logfile 
	echo "with true allele freqs at:" >> $logfile 
	ls -lh $snpsDir/$(basename $snps).g${gen}.r${rep}.round${round}.alleleCts >> $logfile 
	echo "using SNP table:
	$snps
	with ${pctMissing} % missing genotypes" >> $logfile 

	
	if [  $pctMissing  -gt 0 ]; then
	$scriptDir/HAFpipe_wrapper.sh  -t 3,4 -l ${logfile}.${pctMissing}missing \
		-o $HAFDir \
		-s $snps.${pctMissing}pctMissing \
		-m imputed \
		-b ${simreadsStem}.bam \
		-r $refFasta \
		-e $encoding \
		-a $recRate \
		-g $gen \
		-q $quantile
	fi

	if [ "$pctMissing" == 0 ]; then 
	$scriptDir/HAFpipe_wrapper.sh  -t 3,4 -l ${logfile}.${pctMissing}missing \
		-o $HAFDir \
		-s $snps \
		-m none \
		-b ${simreadsStem}.bam \
		-r $refFasta \
		-e $encoding \
		-a $recRate \
		-g $gen	\
		-q $quantile
	fi	


	done
}
export -f calcHAFs

if [ $doCalcHAFs -gt 0 ]; then
	for pctMissing in ${pctMissings[*]}; do

	## add missing sites to SNP table if required
	if [ $pctMissing -gt 0 ] && [ ! -e ${snps}.${pctMissing}pctMissing.numeric.bgz ]; then
		Rscript /mnt/cages/scripts/HAFpipe-line/simulations/addMissingSitesToSNPtable.R $snps $pctMissing
		/mnt/cages/scripts/HAFpipe-line/count_SNPtable.sh ${snps}.${pctMissing}pctMissing
		/mnt/cages/scripts/HAFpipe-line/impute_SNPtable.sh ${snps}.${pctMissing}pctMissing
		/mnt/cages/scripts/HAFpipe-line/prepare_SNPtable_for_HAFcalc.sh ${snps}.${pctMissing}pctMissing
	fi
	echo "run$runID HAFS: parallel --gnu -j$threads calcHAFs ::: ${regimes[*]} ::: ${rounds[*]} ::: ${reps[*]} ::: ${gens[*]} ::: ${covs[*]} ::: $pctMissing ::: $snps ::: $refFasta ::: $recRate" >> $logfile
	parallel --gnu -j$threads calcHAFs ::: ${regimes[*]} ::: ${rounds[*]} ::: ${reps[*]} ::: ${gens[*]} ::: ${covs[*]} ::: $pctMissing ::: $snps ::: $refFasta ::: $recRate
	echo "run$runID HAFS: complete" >> $logfile
	done
fi

#!/bin/bash

###############
## DESCRIPTION: bash script to make a recombined SNPtable from a founder table and a forqs run
###############
#1) specify a forqs run (dir) and a replicate and generation to sample from
#2)randomly pick a set of poolSize(=100) [from 1 to popSize(=1000)] individual indexes to include in pooled sample, sort 
#3)read in forqs pop file lines till you get to an individual on pool list
#4)write each chrom {+ and -} to separate file: (basename $founderSNPs).g${gen}.r${rep}.ind${ind}+/-
#5)parse chrom to breakpoints/chunks
#for each
#	get founder ID, pos
#	tabix snpFile chrom:pos-nextPos|cut $fID >> outFile
#end
#6) when all inds/chroms have been written, then paste them all as separate columns into new recombined-snptable file
#7) keep track of inds and breakpoints used in a brkpts file
################

###############
## FILES
###############
founderSNPFile=${1}
forqsDir=${2}
gen=${3}
rep=${4}
poolSize=${5}
outDir=${6}
chromIX=${7:-1}
forqsPopFile=${forqsDir}/population_gen${gen}_pop${rep}.txt

###############
## SET OUTPUT LOCATIONS
###############
if [ "$outDir" ]; 
then outSNPs=${outDir}/$(basename $founderSNPFile).g${gen}.r${rep}
else outSNPs=/dev/stdout; fi
if [ "$outDir" ]; 
then outBrks=${outDir}/$(basename $founderSNPFile).g${gen}.r${rep}.brkpts
else outBrks=${forqsDir}/$(basename $founderSNPFile).g${gen}.r${rep}.brkpts; fi


## index SNP file if not already done
if [ ! -e ${founderSNPFile}.bgz ]; 
then echo bgzipping and indexing $founderSNPFile; /mnt/cages/ref_data/scripts/SNPtable/index_SNPtable.sh $founderSNPFile; fi


#echo -e 'recombining SNPs in \n' $founderSNPFile '\naccording to breakpoints in\n' $forqsPopFile '\nand writing to\n' $outSNPs'\n and \n'$outBrks

###############
## OPTIONS
###############
nofFounders=$(( $(head -1 $founderSNPFile | tr ',' ' ' | wc -w) - 2 ))

###############
## GET CHROM INFO
###############
forqsChromCount=$(head -2 $forqsPopFile  | tail -1 | cut -f2 -d' ')
snpChrom=$(head -1 $founderSNPFile | cut -f1 -d',')
chromLength=$(tail -1 $founderSNPFile | cut -f1 -d',')


###############
## PICK INDIVIDUALS
###############
forqsPopSize=$(head -1 $forqsPopFile  | cut -f2 -d' ')
indIDList=$(shuf -i 1-$forqsPopSize -n $poolSize | sort -k1g)

###############
## SET UP OUTFILES
###############
cat $founderSNPFile | cut -f1-2 -d',' >  $outSNPs.ind0
echo "founderIX,chromName,fstart,fstop,snpCt" > $outBrks

###############
## GO THROUGH POP FILE
##############
tail -n +3 $forqsPopFile | awk -F ' ' -v indIDList="$indIDList" -v nofFounders="$nofFounders" -v chromName="$snpChrom" -v chromIX="$chromIX" -v chromLength="$chromLength" -v snpFile="$founderSNPFile" -v outSNPs="$outSNPs" -v outBrks="$outBrks" '
	BEGIN {
		indCt=split(indIDList,indIDs," ")
		lookingFor=1
	}
	(NF==0){  ## blank line means NEW IND
		if (currentInd==indIDs[lookingFor]) { # if the previous individual was one you were looking for
			lookingFor++; # update search to next ind on list
		}
		currentInd++; currentChromIX=0; # now update currentInd for lines that come next
	}
	($1=="+"){currentChromIX++;} # start with positive strand
	(NF>0 && currentInd==indIDs[lookingFor] && currentChromIX==chromIX) {  # if this is your ind, and your chrom - proceed!
		## FOUND IND FOR POOL
		indName="ind"lookingFor"forqs"currentInd"_"$1
		printHere=outSNPs".ind"lookingFor $1
		#print indName
		printf indName"\n" > printHere
				
		## GET FIRST CHUNK (start and founder)
		split($3,chunkParse,",")
		fstart=substr(chunkParse[1],2)
		forqsIX=substr(chunkParse[2],1,length(chunkParse[2])-1)
		founderIX=((forqsIX-forqsIX%2)/2)%nofFounders +1 
		
		
		## GO THROUGH REST OF THE CHUNKS
		i=4
		while(fstart<chromLength) {
			if (i==NF)  ## last chunk
			{ fstop=chromLength }
			else { # get chunk end
				split($i,chunkParse,",")
				fstop=substr(chunkParse[1],2)-1
				if (fstop > chromLength) { fstop=chromLength }
			}
			## GET SNPS IN CHUNK FROM FOUNDER
			tabixSNPCommand=("tabix -S1 -p vcf "snpFile".bgz "chromName":"fstart"-"fstop" | cut -f"founderIX+2" -d','")
 
			snpCounter=0
			while((tabixSNPCommand  | getline snpLine ) > 0)
			{
				printf snpLine"\n" >> printHere; snpCounter++
			}
			close(tabixSNPCommand)
			#if (snpCounter < 1) {print tabixSNPCommand; exit 1}
			#print "founder"founderIX","chromName":"fstart"-"fstop","snpCounter" snps"
			print founderIX","chromName","fstart","fstop","snpCounter >> outBrks
			
			fstart=fstop+1
			forqsIX=substr(chunkParse[2],1,length(chunkParse[2])-1) 
			founderIX=((forqsIX-forqsIX%2)/2)%nofFounders +1 
			i++
		}		
	}
'
###############
## PASTE INTO ONE BIG FILE
##############
paste ${outSNPs}.ind* -d',' > ${outSNPs}
rm ${outSNPs}.ind*
#Rscript /mnt/cages/scripts/transposeCSV.R $outFile
#mv ${recSNPsFile}.t $recSNPsFile

#!/bin/bash

###############
## DESCRIPTION: bash script to make a recombined SNPtable from a founder table and a forqs run
###############
#1)get a set of poolSize(=100) <1,popSize(=1000)> individual indexes to include in pooled sample, sort 
#2)read in forqsFile lines till you get to one on list
#3)add new line to outFile
#4)parse to breakpoints/chunks
#for each
#	get founder ID, pos
#	tabix snpFile chrom:pos-nextPos|cut $fID >> outFile
#end
#5)use perl to transpose file from genomesxpositions to positionsxgenomes
################

###############
## FILES
###############
founderSNPFile=${1}
forqsDir=${2}
gens=${3}
rep=${4}
poolSize=${5}
outFile=${6:-${founderSNPFile}.recombined}
forqsPopFile=${forqsDir}/population_gen${gens}_pop${rep}.txt

 
echo -e 'recombining SNPs in \n' $founderSNPFile '\naccording to breakpoints in\n' $forqsPopFile '\nand writing to\n' $outFile

###############
## OPTIONS
###############
nofFounders=$(( $(head -1 $founderSNPFile | tr ',' ' ' | wc -w) - 2 ))

###############
## GET CHROM INFO
###############
forqsChromCount=$(head -2 $forqsPopFile  | tail -1 | cut -f2 -d' ')
chromNames="2L 2R 3L 3R X"
chromNames=$(echo $chromNames | cut -f1-$forqsChromCount -d' ')
snpChrom=$(head -1 $founderSNPFile | cut -f1 -d',')
chromIX=$(echo $chromNames | awk -F ' ' -v snpChrom="$snpChrom" '{chromIX=0; for (i=1;i<=NF;i++){if($i==snpChrom){chromIX=i}};print chromIX}')
chromLength=$(tail -1 $founderSNPFile | cut -f1 -d',')

###############
## PICK INDIVIDUALS
###############
forqsPopSize=$(head -1 $forqsPopFile  | cut -f2 -d' ')
indIDList=$(shuf -i 1-$forqsPopSize -n $poolSize | sort -k1g)

###############
## SET UP OUTFILE
###############
positions=$(cat $founderSNPFile | cut -f1 -d',')
ref=$(cat $founderSNPFile | cut -f2 -d',')
echo $positions | tr ' ' ',' > $outFile
echo $ref | tr ' ' ',' >> $outFile
echo "founderIX,chromName,fstart,fstop,snpCt" > $outFile".breakpts"

###############
## GO THROUGH POP FILE
##############
tail -n +3 $forqsPopFile | awk -F ' ' -v indIDList="$indIDList" -v nofFounders="$nofFounders" -v chromName="$snpChrom" -v chromIX="$chromIX" -v chromLength="$chromLength" -v snpFile="$founderSNPFile" -v outFile="$outFile" ' 
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
		print indName
		printf indName >> outFile
		
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
				printf ","snpLine >> outFile; snpCounter++
			}
			close(tabixSNPCommand)
			#if (snpCounter < 1) {print tabixSNPCommand; exit 1}
			print "founder"founderIX","chromName":"fstart"-"fstop","snpCounter" snps"
			print founderIX","chromName","fstart","fstop","snpCounter >> outFile".breakpts"
			
			if (fstop==chromLength)
			{ printf "\n" >> outFile;}
			fstart=fstop+1
			forqsIX=substr(chunkParse[2],1,length(chunkParse[2])-1) 
			founderIX=((forqsIX-forqsIX%2)/2)%nofFounders +1 
			i++
		}		
	}
'
###############
## IMPORT TO R AND TRANSPOSE
##############
Rscript /mnt/cages/scripts/transposeCSV.R $outFile
mv ${outFile}.t $outFile

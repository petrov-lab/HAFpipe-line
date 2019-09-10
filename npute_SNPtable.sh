#!/bin/bash
	
###############3
## HELP
#################
if [ $# -lt 1 ]
then
    echo "format: npute_SNPtable.sh [snptable] [npute_dir] [winsize=20]"; 
    echo 'One or more variables are undefined. Please try again. 
'
exit
fi

###############
## MAIN
################
snptable=${1}
nputedir=${2}
winsize=${3:-20}
mode=${4:-0}
outfile=${5:-$snptable.npute}


tail -n +2 ${snptable} | cut -f3- -d',' | tr 'N' '?' > ${snptable}.nputeIN
if [ "$mode" = 1 ]; then 
	echo "testing window sizes $winsize"
	python $nputedir/NPUTE.py -m $mode -r $winsize -i ${snptable}.nputeIN -o ${snptable}.nputeTestWins$(echo $winsize | tr ':' '_')
	rm ${snptable}.nputeIN 
else
	echo "running npute with window size $winsize" echo $mode
	python $nputedir/NPUTE.py -m $mode -w $winsize -i ${snptable}.nputeIN -o ${snptable}.nputeOUT
	head -1 $snptable > ${snptable}.npute
	paste <(tail -n +2 ${snptable} | cut -f1-2 -d',') <(cat ${snptable}.nputeOUT | tr 'a-z' 'A-Z' | tr '\t' ',') | tr '\t' ',' >> ${snptable}.npute
	rm ${snptable}.nputeIN ${snptable}.nputeOUT
fi

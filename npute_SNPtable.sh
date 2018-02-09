#!/bin/bash
	
###############3
## HELP
#################
if [ $# -lt 1 ]
then
    echo "format: npute_SNPtable.sh [snptable] [winsize=20]"; 
    echo 'One or more variables are undefined. Please try again. 
'
exit
fi

###############
## MAIN
################
snptable=${1}
winsize=${2:-20}

tail -n +2 ${snptable} | cut -f3- -d',' | tr 'N' '?' > ${snptable}.nputeIN
python software/npute/NPUTE.py -m 0 -w $winsize -i ${snptable}.nputeIN -o ${snptable}.nputeOUT
head -1 $snptable > ${snptable}.npute
paste <(tail -n +2 ${snptable} | cut -f1-2 -d',') <(cat ${snptable}.nputeOUT | tr 'a-z' 'A-Z' | tr '\t' ',') | tr '\t' ',' >> ${snptable}.npute
rm ${snptable}.nputeIN ${snptable}.nputeOUT

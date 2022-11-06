#!/bin/bash


# ==================================================================================================
#      Help
# ==================================================================================================

if [ $# -lt 1 ]
then
    echo "format: npute_SNPtable.sh [snptable] [npute_dir] [winsize=20]";
    echo 'One or more variables are undefined. Please try again.

    **NPUTE python scripts can be downloaded from: http://compgen.unc.edu/wp/?page_id=57
    **the versions included with this repository have been slightly edited to accommodate ambiguous base codes
    **note that python 2 and the numpy package are required for npute'
    exit 1
fi

# ==================================================================================================
#      Main
# ==================================================================================================

snptable=${1}
winsize=${2:-20}
mode=${3:-0}
outfile=${4:-$snptable.npute}
nputedir=${5:-"$(dirname $0)/../NPUTEv1"}

tail -n +2 ${snptable} | cut -f3- -d',' | tr 'N' '?' > ${snptable}.nputeIN
if [ "$mode" = 1 ]; then
    echo "testing window sizes $winsize"
    python $nputedir/NPUTE.py -m $mode -r $winsize -i ${snptable}.nputeIN -o ${snptable}.nputeTestWins$(echo $winsize | tr ':' '_')
    rm ${snptable}.nputeIN
else
    echo "running npute with window size $winsize"
    python $nputedir/NPUTE.py -m $mode -w $winsize -i ${snptable}.nputeIN -o ${snptable}.nputeOUT
    head -1 $snptable > ${snptable}.npute
    paste <(tail -n +2 ${snptable} | cut -f1-2 -d',') <(cat ${snptable}.nputeOUT | tr 'a-z' 'A-Z' | tr '\t' ',') | tr '\t' ',' >> ${snptable}.npute
    rm ${snptable}.nputeIN ${snptable}.nputeOUT
fi

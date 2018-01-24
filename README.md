# HAF-pipe 🏂

The workflow for this tool includes: 🏂

1: install harp (see bitbucket.com/dkessner/harp) <br>
2: install R (tested on versions >= 3.2) and make sure it is in your path

Run the following steps separately for each chromosome: <br>
3: VCF-to-SNPtable.sh <br>
----(input: vcf file, chromosome) <br>
----(optional input: founder subset list)<br>
----(output: table, rows are segregating sites, columns are sequenced founders) <br>
4: imputeMissingCalls.sh<br>
----(input: SNPtable)<br>
----(output:imputed SNPtable)<br>
5: calculateHaplotypeFreqs.sh<br>
----(input: list of bam files, imputed snp table, window size)<br>
----(output: hfreq file for each bam file)<br>
6: calculateAlleleFreqs.sh<br>
----(input:list of hfreq files, snptable (unimputed) )<br>
----(output: af file for each bam file)<br>

# better-HAF-than-AF

The workflow for this tool includes:

1: install harp (see bitbucket.com/dkessner/harp) <br>
Run the following steps separately for each chromosome: <br>
2: VCF-to-SNPtable.sh <br>
----(input: vcf file, chromosome) <br>
----(optional input: founder subset list)<br>
----(output: table, rows are segregating sites, columns are sequenced founders) <br>
3: imputeMissingCalls.sh<br>
----(input: SNPtable)<br>
----(output:imputed SNPtable)<br>
4: calculateHaplotypeFreqs.sh<br>
----(input: list of bam files, imputed snp table, window size)<br>
----(output: hfreq file for each bam file)<br>
5: calculateAlleleFreqs.sh<br>
----(input:list of hfreq files, snptable (unimputed) )<br>
----(output: af file for each bam file)<br>

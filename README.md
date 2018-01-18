# better-HAF-than-AF

The workflow for this tool includes:

1: install harp (see bitbucket.com/dkessner/harp)
Run the following steps separately for each chromosome:
2: VCF-to-SNPtable.sh 
----(input: vcf file, chromosome)  
----(optional input: founder subset list)
----(output: table, rows are segregating sites, columns are sequenced founders) 
3: imputeMissingCalls.sh
----(input: SNPtable)
----(output:imputed SNPtable)
4: calculateHaplotypeFreqs.sh
----(input: list of bam files, imputed snp table, window size)
----(output: hfreq file for each bam file)
5: calculateAlleleFreqs.sh
----(input:list of hfreq files, snptable (unimputed) )
----(output: af file for each bam file)

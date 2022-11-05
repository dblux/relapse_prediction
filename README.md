## TO-DO

- [] Clean up evaluation-figures.ipynb
- [] Clean up sandbox.ipynb

## Scripts
1. [Python script](bin/calcGeneLength.py) for calculating mean, median and max gene lengths from reference genome GTF file. Saves output file as gene_lengths.tsv by default

To run:

```	
# To calculate gene lengths from NCBI reference genome GTF
$ python calcGeneLength.py ncbi ref_genome.gtf

# To specify outfile path
$ python calcGeneLength.py ensembl ref_genome.gtf -o output_file.tsv

# To show help message
$ python calcGeneLength.py -h
```

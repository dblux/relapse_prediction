#!/usr/bin/env python
# coding: utf-8

import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input", # Positional argument
                    help = "Ensembl GTF file path")
parser.add_argument("--output", # Optional argument
                    help = "Output file path")
args = parser.parse_args()
print(args.input)

RPATH = "../../info/ref_genome/GRCh38/ensembl-Homo_sapiens.GRCh38.99.gtf"
WPATH = "gene_lengths.tsv"

#%%

class Gene:
    def __init__(self, name):
        self.name = name
        self.transcripts = {}
    
    def calcMean(self):
        transcript_lengths = [
            sum(list_exons) for list_exons in self.transcripts.values()
        ]
        return np.mean(transcript_lengths)
    
    def calcMedian(self):
        transcript_lengths = [
            sum(list_exons) for list_exons in self.transcripts.values()
        ]
        return np.median(transcript_lengths)
    
    def calcMax(self):
        transcript_lengths = [
            sum(list_exons) for list_exons in self.transcripts.values()
        ]
        return max(transcript_lengths)

#%%
def main():
    fw = open(WPATH, "w")
    header = "gene_name\tmean_length\tmedian_length\tmax_length\n"
    fw.write(header)
    
    with open(RPATH) as gtf:
        # Assumption: Feature type is always ordered "gene, transcript, exon"
        for line in gtf: # iterator
            if not line.startswith("#"): # ignore comments
                record = line.strip().split("\t")
                
                if record[2] == "gene":
                    # If object gene exists -> Write lengths to file
                    try:
                        print(gene.transcripts)
                        wline = "{}\t{}\t{}\t{}\n".format(
                            gene.name, gene.calcMean(),
                            gene.calcMedian(), gene.calcMax()
                        )
                        fw.write(wline)
                        # break
                    except NameError: # catch first error: no gene object
                        pass
                    
                    # Instantiate gene
                    gene_name = record[8].split(";")[2].split('"')[1]
                    gene = Gene(gene_name)
                    
                elif record[2] == "transcript":
                    # Add (transcript name: list) to dict
                    transcript_name = record[8].split(";")[7].split('"')[1]
                    assert transcript_name.startswith(gene.name)
                    gene.transcripts[transcript_name] = []
                    
                elif record[2] == "exon":
                    exon_gene = record[8].split(";")[5].split('"')[1]
                    assert exon_gene == gene.name
                    # Append length of exon to list
                    exon_length = int(record[4]) - int(record[3]) + 1
                    # Depends on previous transcript name
                    gene.transcripts[transcript_name].append(exon_length)
    
    fw.close()
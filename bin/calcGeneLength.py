#!/usr/bin/env python
# coding: utf-8

import numpy as np
import argparse, re

class Gene:
    def __init__(self, name):
        self.name = name
        self.transcripts = {} # gene_name -> list_exon
        self.synonyms = []
    
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

# Genes with same gene name but different gene IDs are combined
# Return dictionary of objects
def readEnsembl(rpath):
    dict_genes = {} # Initialise dict (gene_name -> gene_obj)
    
    def save2Dict(gene):
        if gene.name not in dict_genes.keys():
            # If transcripts dict is not empty add gene object to dict
            if len(gene.transcripts) > 0:
                dict_genes[gene.name] = gene
            
        else:
            # Check that keys in both dict do not intersect
            # before merging
            assert not any(
                key in dict_genes[gene.name].transcripts
                for key in gene.transcripts
            )
            
            dict_genes[gene.name].transcripts.update(
                gene.transcripts
            )
    
    gtf = open(rpath, "r")
    # Assumption: Feature type is always ordered "gene, transcript, exon"
    for line in gtf: # iterator
        if not line.startswith("#"): # ignore comments
            record = line.strip().split("\t")
            
            if record[2] == "gene":
                try: # Try to save gene
                    save2Dict(gene)
                        
                except NameError: # catch first error: no gene object
                    pass
                
                # Instantiate new gene object
                gene_name = record[8].split(";")[2].split('"')[1]
                gene = Gene(gene_name)
                
            elif record[2] == "transcript":
                # Add (transcript name: list) to dict
                transcript_name = record[8].split(";")[7].split('"')[1]
                assert transcript_name.startswith(gene.name)
                # Check that key is not already present in dict
                assert transcript_name not in gene.transcripts.keys()
                gene.transcripts[transcript_name] = [] # initialise empty list
                
            elif record[2] == "exon":
                exon_transcript = record[8].split(";")[8].split('"')[1]
                assert exon_transcript == transcript_name
                # Append length of exon to list
                exon_length = int(record[4]) - int(record[3]) + 1
                # Depends on previous transcript name
                gene.transcripts[transcript_name].append(exon_length)
    
    # Saves last gene object
    save2Dict(gene)
    gtf.close() 
    
    return dict_genes

# Genes with same gene name but different gene IDs are combined
# Return dictionary of objects
def readNCBI(rpath):
    dict_genes = {} # Initialise dict (gene_name -> gene_obj)
    previous_transcript = "dummy"
    
    def save2Dict(gene):
        if gene.name not in dict_genes.keys():
            # If transcripts dict is not empty add gene object to dict
            if len(gene.transcripts) > 0:
                dict_genes[gene.name] = gene
            
        else:
            # Check that keys in both dict do not intersect
            # before merging
            assert not any(
                key in dict_genes[gene.name].transcripts
                for key in gene.transcripts
            )
            
            dict_genes[gene.name].transcripts.update(
                gene.transcripts
            )
    
    gtf = open(rpath, "r")
    # Assumption: Feature type is always ordered "gene, transcript, exon"
    for line in gtf: # iterator
        if not line.startswith("#"): # ignore comments
            record = line.strip().split("\t")
            
            if record[2] == "gene":   
                # Save previous gene object to dict
                try: # try to save gene
                    save2Dict(gene)
                    
                except NameError: # catch first error: no gene object
                    pass
                
                # Instantiate new gene object
                gene_name = record[8].split(";")[0].split('"')[1]
                gene = Gene(gene_name)
                # Add synonyms
                gene.synonyms = re.findall("gene_synonym \"(.*?)\"", record[8])
                
                if record[1] == "Curated Genomic":
                    # Add "exon" length
                    exon_length = int(record[4]) - int(record[3]) + 1
                    gene.transcripts["curated_genomic"] = [exon_length]
                
            elif record[2] == "exon":
                exon_transcript = record[8].split(";")[1].split('"')[1]
                exon_length = int(record[4]) - int(record[3]) + 1
                
                if exon_transcript != previous_transcript:
                    # Check that transcript does not appear >1 in each gene
                    assert exon_transcript not in gene.transcripts.keys()
                    # Create new list and assign exon length
                    gene.transcripts[exon_transcript] = [exon_length]
                    previous_transcript =  exon_transcript
                    
                else:
                    # Append length of exon to list    
                    gene.transcripts[exon_transcript].append(exon_length)
    
    # Saves last gene object
    save2Dict(gene)
    gtf.close()
    
    return dict_genes

# Write lengths to file
def writeEnsembl(dict_genes, wpath):
    f = open(wpath, "w")
    header = "gene_name\tmean_length\tmedian_length\tmax_length\n"
    f.write(header)
    
    for gene in dict_genes.values():        
        wline = "{}\t{}\t{}\t{}\n".format(
            gene.name, gene.calcMean(),
            gene.calcMedian(), gene.calcMax()
        )
        f.write(wline)
    
    f.close()

# Write lengths to file
def writeNCBI(dict_genes, wpath):
    f = open(wpath, "w")
    header = "gene_name\tmean_length\tmedian_length\tmax_length\tgene_synonyms\n"
    f.write(header)
    
    for gene in dict_genes.values():        
        list_synonyms = ",".join(gene.synonyms)
        wline = "{}\t{}\t{}\t{}\t{}\n".format(
            gene.name, gene.calcMean(),
            gene.calcMedian(), gene.calcMax(), list_synonyms
        )
        
        f.write(wline)
        
    f.close()
   
if __name__ == "__main__":
    DOCSTRING = "Calculates gene lengths from reference genome GTF file"
    WPATH = "gene_lengths.tsv"
    
    parser = argparse.ArgumentParser(description=DOCSTRING)
    parser.add_argument("format",
                        choices = ["ncbi", "ensembl"],
                        help = "GTF file format")
    parser.add_argument("infile",
                        help = "Input GTF file path")
    parser.add_argument("-o", "--outfile",
                        default=WPATH,
                        help = "Output file path")
                        
    args = parser.parse_args()
    
    if args.format == "ncbi":    
        dict_genes = readNCBI(args.infile)
        writeNCBI(dict_genes, args.outfile)
        
    elif args.format == "ensembl":
        dict_genes = readEnsembl(args.infile)
        writeEnsembl(dict_genes, args.outfile)

    # for obj in dict_genes.values():
    #     print(obj.name, obj.transcripts)

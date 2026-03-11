# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 15:41:13 2026

@author: nilof
"""
"""
This script uses gtf files, extracts the relevent information and 
stores it in the genome and gene classes
"""

from genome_gene_classes import Genome, Gene


def extract_genome_info(genome_file):
    genome_text = open(genome_file, "r")  #reading the file
    genome_id=""
    
    
    #get the genome file ID
    for i, line in enumerate(genome_text, start= 1):
        if i==3:
            genome_id= line
            break
    
    genome= Genome(genome_id)
    print("Genome object created:", genome.genome_id)
    current_gene= None
    
    for line in genome_text:
        fields= line.strip().split("\t")
        if len(fields) < 3:    #skip the intro line
            continue 
        
        feature= fields[2]
        start= int(fields[3])
        end= int(fields[4])
        
        # --- CASE 1: new gene starts ---
        if feature == "gene":
            # save previous gene
            if current_gene is not None:
                genome.add_gene(current_gene)
                

            # create new gene
            gene_name = fields[8]  # or parse ID=...
            current_gene = Gene(gene_name, start, end)
            
        # --- CASE 2: exon inside current gene ---  # 
        elif feature == "exon" and current_gene is not None:
            current_gene.exons.append((start, end))

    # save last gene
    if current_gene is not None:
        genome.add_gene(current_gene)
        
    print("number of gene objects created", genome.get_num_genes())
    return genome
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 15:41:16 2026

@author: nilof
"""
### a script to test classes and functions imported below
### uses Guilt genome, which is the smallest genome in eaukaryotes 

from extract import extract_genome_info
from genome_gene_classes import Genome, Gene
import matplotlib.pyplot as plt

A= extract_genome_info("genomic.gtf")


#creates a histogram with number of giant introns of a genome
#returns the giant_introns_list
def get_genome_giant_introns(genome_object, threshold):
    giant_introns = []
    for gene in genome_object.genes.values():
        
        #some lines for debugging
        # print("gene_start and gene_end", gene.start, gene.end)
        # print(gene.exons)
        
        gene.calculate_intron_lenghts(gene.exons)
        intron_list = gene.get_intron_lenghts()
        # print(intron_list) #for debugging
        for item in intron_list:
            if item > threshold:
                giant_introns.append(item)
                
    plt.hist(giant_introns, int(max(giant_introns)/1000))           
    plt.show()
    return giant_introns
    
    
          
B= get_genome_giant_introns(A, 3000)
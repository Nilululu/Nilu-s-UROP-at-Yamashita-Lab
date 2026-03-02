# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 16:00:10 2026

@author: nilof
"""

### We will use the following classes to rpoperly store relevent data

class Gene:
    
    """
    CLASS GENE stores a gene name, position, introns and exons
    """
    def __init__(self, name, start, end, exons= None):
        
        self.name = name
        self.start = start
        self.end = end
        self.length = end-start
        self.exons= exons.sort(key=lambda x: x[0])if exons else []
        self.intron_lengths= []
        
        
    #calculating intron lenght from exons positions    
    def calculate_intron_lenghts(self, exons):
        if len(exons)<= 1:
            return 
        
        else:
            exons.sort(key=lambda x: x[0])
            for i in range(len(exons)-1):
                intron_start= int(exons[i][1])
                intron_end= int(exons[i+1][0])
                intron_lenght= intron_end - intron_start
                self.intron_lengths.append(intron_lenght)
            
            # return 
    
    #fuctions to retrive intron and exons information
    def get_intron_lenghts(self):
        return self.intron_lengths
    
    def get_sorted_intron_lens(self):
        intron_list= list(self.intron_lenght)
        return intron_list.sort()



class Genome:
    
    """
    CLASS GENOME creates a genome object, can add and retireve genens
    """
    def __init__(self, genome_id, genes= None):
        self.genome_id= genome_id
        self.genes= dict (genes) if genes else {}
        
    def add_gene (self, gene):
        self.genes[gene.name]= gene
        
    def get_gene (self, name):
        return self.genes[name]
    
    def get_all_genes (self):
        return self.genes
    
    def get_num_genes (self):
        return len(self.genes)
    
    


    
    
    
    
    
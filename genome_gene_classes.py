# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 16:00:10 2026

@author: nilof
"""

class Gene:
    def __init__(self, name, start, end, axons= None):
        
        self.name = name
        self.start = start
        self.end = end
        self.length = end-start
        self.axons= list(axons) if axons else []
        self.intron_lengths= []
        
        
        
    def calculate_intron_lenghts(self, axons):
        if len(axons)< 1:
            return 
        
        else:
            for i in range(len(axons)-2):
                intron_start= int(axons[i][1])
                intron_end= int(axons[i+1][0])
                intron_lenght= intron_end - intron_start
                self.intron_lengths.add(intron_lenght)
            
            return 
        
    def get_intron_lenghts(self):
        return self.intron_lengths
    
    def get_sorted_intron_lens(self):
        intron_list= list(self.intron_lenght)
        return intron_list.sort()
    
class Genome:
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
    
    



# def genome_info(genome_file):
#     genome_text = open(genome_file, "r") 
#     genome_id=""
    
    
#     #get the genome file ID
#     for i, line in enumerate(genome_text, start= 1):
#         if i==3:
#             genome_id= line
#             break
    
#     genome= Genome(genome_id)
#     current_gene= None
    
#     for line in genome_text:
#         fields= line.strip().split("\t")
#         if len(fields) < 3:    #skip the intro line
#             continue 
        
#         feature= fields[2]
#         start= int(fields[3])
#         end= int(fields[4])
        
#         # --- CASE 1: new gene starts ---
#         if feature == "gene":
#             # save previous gene
#             if current_gene is not None:
#                 genome.add_gene(current_gene)

#             # create new gene
#             gene_name = fields[8]  # or parse ID=...
#             current_gene = Gene(gene_name, start, end)

#         # --- CASE 2: exon inside current gene ---  # 
#         elif feature == "exon" and current_gene is not None:
#             current_gene.axons.add((start, end))

#     # save last gene
#     if current_gene is not None:
#         genome.add_gene(current_gene)

    
    
    
    
    
    
    
    
    
    
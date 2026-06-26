#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 11:08:23 2026 

Scirpt for getting general genomic folder data and writing it in results_table.txt, contains 
    introns stats: extracts speciefic statistics about introns and return them in for of a list
    write_to_table: exctracts the data of interest from genomic folder and returns it in form of a string
    multiprocessor used to loop through all genomic files using write to table
    Read write_to_table to learn about what data is being stored in the text file created
"""

#pyhton modules 
import logging
from multiprocessing import Pool
import numpy as np
from pathlib import Path

#project modules 
from extract import extract_id_and_genes, compute_intron 
from metadata import get_genome_metadata
import taxonomy 


logging.basicConfig(filename="info_test1_logger.txt", level=logging.DEBUG, force=True) 
logger = logging.getLogger(__name__)
logger.info("the loggin initialized successfully")



# creating a taxonomy dictionary from taxonomy database

nodes_file = "nodes.dmp"
names_file = "names.dmp"
taxonomy_dict = taxonomy.generate_taxonomy_dict(nodes_file)
tax_to_name = taxonomy.generate_tax_to_name(names_file)




def intron_stats (genes):
    """
    helper function for write_to_table, calculates statistical values about 
    introns from a genome dictionary 

    Parameters
    ----------
    genes : a dictionary (already includes intron start and end length insides)

    Returns : 
    -------
    max_intron, min_intron, mean_intron,median_intron, sd_intron, q_25, q_50, 
    q_75, q_95, q_99, q_999, q_9999, q_99999
    """
    introns = []
    
    for gene in genes:
        for item in genes[gene]['introns']:
            if item != set():  #skipping genes with no introns
                start = item[0]
                end = item [1]
                length = end - start
                # if length < 0 :
                #     print(length)
                                            
                introns.append(length)
    
    if introns: 
        max_intron = max(introns)
        min_intron = min(introns)
        mean_intron = int(np.mean(introns))
        median_intron = int(np.median(introns))
        sd_intron = int(np.std(introns))
        
        q_25 = int(np.quantile(introns, 0.25))
        q_50 = int(np.quantile(introns, 0.5))
        q_75 = int(np.quantile(introns, 0.75))
        q_95 = int(np.quantile(introns, 0.95))
        q_99 = int(np.quantile(introns, 0.99))
        q_999 = int(np.quantile(introns, 0.999))
        q_9999 = int(np.quantile(introns, 0.9999))
        q_99999  = int(np.quantile(introns, 0.99999))       
    
        intron_stats = [max_intron, 
            min_intron, mean_intron,median_intron, sd_intron, q_25, q_50, q_75, q_95, q_99, q_999, q_9999, q_99999
            ]

    else:
        intron_stats = [0,0,0,0,0,0,0,0,0,0,0,0,0]
        
    return intron_stats








    
def write_to_table (line):
    """
    exctracts the data of interest from a gtf file and returns it in the following order as a tab saperated string
    
    genome_id, name, kingdom, tax_id, total_sequence_length, assembly_level, assembly_type, 
    numChr, num_scaffolds, num_contigs, scaffold_n50, contig_n50, gc_percent, max_intron, 
    min_intron, mean_intron,median_intron, sd_intron, q_25, q_50, 
    q_75, q_95, q_99, q_999, q_9999, q_99999
        
    
    Parameters
    ----------
    line : a line from genomic_directory.csv containing the gtf file location

    Returns
    -------
    table_str: str
    """
    gtf_file= Path(line.split(",")[0])
    genome_id, genes = extract_id_and_genes(gtf_file)
    compute_intron(genes)
    
    metadata = get_genome_metadata(gtf_file)  #information about genome assembly inlcuding tax num,  sequence length, etc
    
    tax_dict = taxonomy.find_taxonomy(metadata[0], taxonomy_dict, tax_to_name)
    kingdom = tax_dict["kingdom"]
    name = tax_dict["species"]
    
    
    intron_stats_list = intron_stats(genes)
        
    table_list = [genome_id, name, kingdom]
    table_list.extend(metadata)
    table_list.extend(intron_stats_list)
    
    
    table_str = ""
    
    for element in table_list:
        table_str = table_str + str(element) + "\t"
    
    return table_str
        

    

genomic_directory = "genomic_directory.csv"

with open(genomic_directory, 'r') as directory:
    lines = directory.readlines()
    
    if __name__ == '__main__':

        with Pool(5) as p:
            
            
            results = p.map(write_to_table, lines)
     
            with open("result_table.txt", 'w') as table:
                
                # adding column names to the table
                info_line = "#genome_id name kingdom tax_id total_sequence_length assembly_level assembly_type numChr num_scaffolds num_contigs scaffold_n50 contig_n50 gc_percent max_intron min_intron mean_intron median_intron sd_intron q_25 q_50 q_75 q_95 q_99 q_999 q_9999 q_99999"
                info_line.replace(" ", "\t")
                table.write(info_line)
                table.write("\n")
                
                #writing the data 
                for item in results:
                    table.write(item)
                    table.write("\n")
           




    
    
    
    
    
    
    
    
    

    
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 11:08:23 2026 

Scirpt for parralal analysis of genomic folder data, composed of two main function, 
    write_to_table: exctracts the data of interest from genomic folder and writes it in a txt file 
    mass_reading: loops over a group of genomic folders with write _to_table function 
    

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








    
def write_to_table (gtf_file, table):
    """
    exctracts the data of interest from a genomic folder and writes it in a txt file 
    Writes following information about each genome in a line
    
    genome_id, name, kingdom, tax_id, total_sequence_length, assembly_level, assembly_type, 
    numChr, num_scaffolds, num_contigs, scaffold_n50, contig_n50, gc_percent, max_intron, 
    min_intron, mean_intron,median_intron, sd_intron, q_25, q_50, 
    q_75, q_95, q_99, q_999, q_9999, q_99999
        
    
    Parameters
    ----------
    gtf_file : srt file path to genomic.gtf
    table : txt file for writing the data

    Returns
    -------
    None
    """
    
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
    
    with open (table, 'a') as open_table:
        open_table.write(table_str + "\n")
    return
        
       
        
       
        
       
        
       
        
def mass_reading (n):
    """
    loops over a group of genomic files with write _to_table function 

    Parameters
    ----------
    genomic_directory : a csv file containing the location of all gtf files
    n : list of two integer, start, end indicies to read from in genomic_directory  

    Returns
    -------
    None
    """
    
    genomic_directory = "genomic_directory.csv"
    with open(genomic_directory, 'r') as directory:
        lines = directory.readlines()
        
        start = n[0]
        end = n[1]
        table = "table_{}_{}.txt".format(str(start), str(end))
        
        for line in lines[start:end]:
            try:
                
                gtf_loc= Path(line.split(",")[0])
                write_to_table(gtf_loc, table)
                
            except:
                logger.error(f"failed to extract information for {line}")
                raise
    logger.info(table)      
    return





index_list = list(range(0,2641, 264))  # for cluster

# index_list = list(range(0, 31, 6))   # for testing on my laptop

for i in range(len(index_list)-1):
    index_list[i] = [index_list[i], index_list[i+1]]

index_list = index_list[:-1]    
    

    
if __name__ == '__main__':

    with Pool(10) as p:
        p.map(mass_reading, index_list)


    
    
    
    
    
    
    
    
    

    
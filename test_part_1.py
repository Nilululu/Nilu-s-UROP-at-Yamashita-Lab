#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 11:08:23 2026

@author: nilofarghafory
"""
import logging



logging.basicConfig(filename="info_test1_logger.txt", level=logging.DEBUG, force=True) 
logger = logging.getLogger(__name__)
logger.info("the loggin initialized successfully")

from extract import extract_genesAndId, compute_intron 
from metadata import get_GintronInfo, get_genomeMetadata
from pathlib import Path
# from plot_creator import create_hist, create_scatter, create_2d_scatter, get_style
# import numpy as np
import taxonomy 
# import matplotlib.pyplot as plt
# import math


nodes_file = "nodes.dmp"
names_file = "names.dmp"

taxonomy_dict = taxonomy.generate_taxonomy_dict(nodes_file)
tax_to_name = taxonomy.generate_tax_to_name(names_file)


def write_to_table (gtf_file, table, gintron_treshold):
    #so for every plot I need 1-3 values, I can start by making a function for one plot and see how to generalize it to everything
    #I can have some lists, I just need to empty them after each use, or each time I call the function
    #I might collect taxNum of all genomes and find their taxonomy all togehter at once 
    # I need to make a table with all the values I need, I can make the table tab saperated and the elements within two tabs be comma saperated?
    #let me draft what kind of table it will be
    
    ###genes_dic and taxonom dictionary are not here, but they will be created and used in the process 
    ###genome_id name status type size numchr taxnum max_intron kingdom
    
    genomeId, genes = extract_genesAndId(gtf_file)
    compute_intron(genes)
    
    Gintron, GintronStart, GintronEnd = get_GintronInfo(genes, genomeId, gintron_treshold)
    
    taxNum, genome_size, status, numChr, genome_type= get_genomeMetadata(gtf_file)
    taxId = taxonomy.find_taxonomy(taxNum, taxonomy_dict, tax_to_name)
    kingdom = taxId["kingdom"]
    name = taxId["species"]
    
    
    max_intron = 0
    
    if Gintron:
        max_intron = max(Gintron)
    
    else:
        for gene in genes:
            for item in genes[gene]['introns']:
                if item != set():  #skipping genes with no introns
                    start = item[0]
                    end = item [1]
                    length = end - start
                                                
                    if length > max_intron:
                        max_intron = length
                    
                else:
                    continue
        
    
    table_list = [genomeId, name, status, genome_type, genome_size, numChr, taxNum, max_intron, kingdom]
    table_str = ""
    
    for element in table_list:
        table_str = table_str + str(element) + "\t"
    
    print(table_str)
    
    with open (table, 'a') as open_table:
        open_table.write(table_str + "\n")
    return
        
        
def mass_reading (genomic_directory, table):
    
    with open(genomic_directory, 'r') as directory:
        lines = directory.readlines()
        for line in lines:
            try:
                
                gtf_loc= Path(line.split(",")[0])
                write_to_table(gtf_loc, table, 1000)
                
            except:
                logger.error(f"failed to extract information for {line}")
                raise
            
    return

mass_reading("genomic_directory.csv", "table.txt")




        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
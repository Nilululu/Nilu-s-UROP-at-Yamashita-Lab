
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 11:08:23 2026 

Script for getting all introns in our genome folders 
    The main called multiprocessor on get_id_and_introns.

    get_id_and_introns: extracts genome_id and genes from a genomic file, and returns the id and all the introns of the genome in a list.
    Multiprocessor is used to loop through all genomic files using get_id_and_introns and record the output in introns.txt file

"""

#pyhton modules 
import logging
from multiprocessing import Pool
from pathlib import Path
import time

#project modules 
from extract import extract_id_and_genes, compute_intron 

### making a logging file
logging.basicConfig(filename="info_test3_logger.txt", level=logging.INFO) 
logger = logging.getLogger(__name__)


def get_id_and_introns (line):
    """
    exctracts all the intron lenghts (accounting for multiple transcrips)from a gtf file and 
    writes it in a txt file
    
    Parameters
    ----------
    line : a line from genomic_directory.csv containing the gtf file location
    
    Returns
    -------
    None
    """
    gtf_file= Path(line.split(",")[0])
    genome_id, genes = extract_id_and_genes(gtf_file)
    compute_intron(genes)
    
    id_and_introns = [genome_id]   # to keep track of genome id and introns 
    
    #calculating intron lengths 
    for gene in genes:
        for item in genes[gene]['introns']:
            if item != set():  #skipping genes with no introns
                start = item[0]
                end = item [1]
                length = end - start
                
                id_and_introns.append(length)
    
    return id_and_introns
    

#using multiprocessor 
    
if __name__ == '__main__':
    
    #getting all the lines for workers for loop through
    genomic_directory = "genomic_directory.csv"

    with open(genomic_directory, 'r') as directory:
        lines = directory.readlines()
    
    start_time = time.time()

    #passing the lines
    with Pool(5) as p:
        results = p.map(get_id_and_introns, lines)
    

        
    #recording the output
    with open("introns.txt", "w") as introns_file:
        
        for info in results:
            info = list(map(str, info))
            info = ("\t").join(info)
            introns_file.write("{}\n".format(info))
    
    end_time = time.time()
    duration = end_time - start_time
    
    #recording the time and duration
    logger.info("start time: {}".format(start_time))   
    logger.info("end time: {}".format(end_time))
    logger.info("duration: {}".format(duration))
    
            
               
                
                
        
        
        
 


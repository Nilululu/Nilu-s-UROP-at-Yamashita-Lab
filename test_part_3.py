
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 11:08:23 2026 

Scirpt for getting all introns in our genome folders 
    write_to_table: exctracts the data of interest from genomic folder and writes it in text files name wokrer_{worker_id}.txt file
    multiprocessor used to loop through all genomic files using write to table

"""

#pyhton modules 
import logging
from multiprocessing import Pool
from pathlib import Path
import os

#project modules 
from extract import extract_id_and_genes, compute_intron 

### making a logging file
logging.basicConfig(filename="info_test3_logger.txt", level=logging.DEBUG) 
logger = logging.getLogger(__name__)
worker_ids = set()


def write_to_table (line):
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
    
    introns = []   # to keep tracck of introns 
    
    #calculating intron lengths 
    for gene in genes:
        for item in genes[gene]['introns']:
            if item != set():  #skipping genes with no introns
                start = item[0]
                end = item [1]
                length = end - start
                # if length < 0 :      #for debugging 
                #     print(length)
                                            
                introns.append(length)
    
    #writing genome id followed by all intron lengths for a genome in one line
   
    worker_id = os.getpid()
    worker_ids.add(worker_id)
    text_file = "worker" + str(worker_id) + ".txt"
    with open(text_file, "a") as table:   #how to change this so they all don't write in the same text file
        table.write(genome_id)
        table.write(" ")
        table.write(str(introns).strip("[").strip("]").replace(", ", " "))
        table.write("\n")
        
    

genomic_directory = "genomic_directory.csv"

with open(genomic_directory, 'r') as directory:
    lines = directory.readlines()
    
    if __name__ == '__main__':

        with Pool(5) as p:
            p.map(write_to_table, lines)
   
logger.info(worker_ids)



# ### plots with all introns 
# ### I want to see how will this violin plot look like for different kingdoms 
# with open("intron_table.txt", 'r') as intron_table:
#     all_introns = []

#     for line in intron_table:
        
#         line_introns = line.split(",")
#         line_introns = [int(x) for x in line_introns]
#         all_introns.extend(line_introns)
    
    
#     introns = np.array(all_introns)
#     fig9, ax9 = plt.subplots()
#     ax9.violinplot(np.log10(introns), showmedians = True, showmeans = True)
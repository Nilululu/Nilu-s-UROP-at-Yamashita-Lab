#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

A script to go through multiple gtf files and exctract information about their genes and introns

Information it logs about each gene
gene_id genome_id kingdom phylum species gene_len non_coding_len ratio max_intron max_intron_start max_intron_end tax_id total_sequence_length assembly_status assembly_level assembly_type numChr num_scaffolds num_contigs scaffold_n50 contig_n50 gc_percent"

Still under work ....

Created on Fri Jul 17 16:34:42 2026 by nilu
"""


#pyhton modules 
import logging
from multiprocessing import Pool
from pathlib import Path
import time

#project modules 
from extract import extract_id_and_genes, compute_intron 
from metadata import get_genome_metadata
import taxonomy 


logging.basicConfig(filename="info_test4_logger.txt", level=logging.INFO, force=True) 
logger = logging.getLogger(__name__)
logger.info("the loggin initialized successfully")



# creating a taxonomy dictionary from taxonomy database

nodes_file = "nodes.dmp"
names_file = "names.dmp"
taxonomy_dict = taxonomy.generate_taxonomy_dict(nodes_file)
tax_to_name = taxonomy.generate_tax_to_name(names_file)





# a function to go in gene and find the gene length, non_coding_sequence_lenght (introns length), 
#ratio of introns length to gene lengths, as well as the biggest intron in a gene and it's relative 
# normalized position with respect to the gene

def gene_info(gene):
    gene_start, gene_end = gene['position']
    gene_len = gene_end - gene_start
    
    non_coding = 0
    max_intron_gene = 0
    max_start, max_end = 0, 0    #relative and normalized to gene position and len
    
    for intron in gene["intron"]:
        if intron:
            start, end = intron
            intron_len = end - start
            non_coding = non_coding + intron_len
            if intron_len > max_intron_gene:
                max_intron_gene = intron_len
                max_start, max_end = (start - gene_start)/gene_len, (end - gene_start)/gene_len
    
    ratio = non_coding/gene_len
    
    output = [gene_len, non_coding, ratio, max_intron_gene, max_start, max_end]
    return output


#now I need a function that takes a genome and go through all the genes in the genome,
# and writes the information in a text file



    
def write_to_table (line):
    
    try: 
        gtf_file= Path(line.split(",")[0])
        genome_id, genes = extract_id_and_genes(gtf_file)
        compute_intron(genes)
        
        metadata = get_genome_metadata(gtf_file)  #information about genome assembly inlcuding tax num,  sequence length, etc
        tax_dict = taxonomy.find_taxonomy(metadata[0], taxonomy_dict, tax_to_name)
        kingdom = tax_dict["kingdom"]
        phylum = tax_dict.get("phylum", "NA")
        species = tax_dict["species"]
        
        all_genes = []
        
        for gene in genes:
            gene_id = gene
            info = gene_info(genes[gene])
            
            gene_data = [gene_id, genome_id, kingdom, phylum, species]
            gene_data.extend(info)
            gene_data.extend(metadata)
            
            all_genes.append(gene_data)
            
        
        return all_genes
    except:
        logger.error("failed to extract information for the following genome {}".format(line))


start_time = time.time()  

genomic_directory = "genomic_directory.csv"  #for refseq annotated ones
#genomic_directory = "genomic_directory_gca.csv"  #for independant annotations 

with open(genomic_directory, 'r') as directory:
    lines = directory.readlines()
    
    if __name__ == '__main__':

        with Pool(5) as p:
            
            #using 5 walkers to go over all the genomes in genomic_direcctory lines
            results = p.map(write_to_table, lines)
            
            
        #writing the results in a text file
        with open("gene_table_refseq.txt", 'w') as table:   # !!!change to gene_table_gca for independent annotations
            
            # adding column names to the table
            info_line = "#gene_id genome_id kingdom phylum species gene_len non_coding_len ratio max_intron max_intron_start max_intron_end tax_id total_sequence_length assembly_status assembly_level assembly_type numChr num_scaffolds num_contigs scaffold_n50 contig_n50 gc_percent"
            info_line = info_line.replace(" ", "\t")
            table.write("{}\n".format(info_line))

            #writing the data 
            for genome in results:
                if genome: 
                    for gene in genome:
                        table.write("{}\n".format(gene))
           

        end_time = time.time()
        duration = end_time - start_time
        
        logger.info("start time: {}".format(start_time))
        logger.info("end time: {}".format(end_time))
        logger.info("duration: {}".format(duration))


    
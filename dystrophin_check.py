#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 14:23:23 2026

This script will go through all NCBI refeq genomes that we have and filter for mammals. 
Afterwards, it will try to find DMD (dystrophin) gene in each mammals and record gene information
The script takes note of mammals with no DMD gene as well. 
The findings are output in a table named "mammals_dystrophin.txt"

"""

import logging


logging.basicConfig(filename= "log_dystrophin.txt", level = logging.ERROR, force = True)
logger = logging.getLogger(__name__)

#importing internal modules
from metadata import get_genome_metadata
from pathlib import Path
from taxonomy import generate_tax_to_name, generate_taxonomy_dict, find_taxonomy
from extract import extract_id_and_genes, compute_intron


#generate taxonomic dictionary
names = "names.dmp"
nodes = "nodes.dmp"

tax_to_name = generate_tax_to_name(names)
taxonomy_dict = generate_taxonomy_dict(nodes)


# using genomic dirsctory 
with open ("genomic_directory.csv", 'r') as directory, open("mammals_dystrophin.txt", 'w') as table:
    
    header = "genome_id tax_id species kingdom db_xref gene_length products max_intron"
    header = ("\t").join(header.split(" ")) 
       
    table.write("{}\n".format(header))
    for line in directory:
        fields = line.split(",")
        loc = Path(fields[0])
        
        #filtering for mammals
        meta = get_genome_metadata(loc)
        status = meta[2]
        if status != "current":
            logger.error("suppressed genome found: {}".format(loc))
            continue
                      
        tax_id = meta[0]
        taxa = find_taxonomy(tax_id, taxonomy_dict, tax_to_name, {})

        taxa_class = taxa.get("class", "No Class")
        if taxa_class == "mammals": 
            #extracting all information for the genome
            genome_id, genome = extract_id_and_genes(loc)
            compute_intron(genome)
            
            DMD_found = False
            
            #searching for DMD
            for gene in genome:
            
                if gene.upper() == "DMD":
                    if DMD_found: #there shouldn't be two dystrophin gene on the same genome
                        logger.error("two DMD gene found on the same genome, {}".format(genome_id))
                    
                    DMD_found = True
                    
                    start, end = genome[gene]["position"]
                    gene_length = end - start
                    
                    products = set()
                    for product in genome[gene]["products"]:
                        if product.startswith("dystrophin"):
                            products.add("dystrophin")
                        else:
                            products.add(product)

                    introns = genome[gene]["introns"]
                    
                    intron_lens = []
                    for intron in introns:
                        istart, iend = intron
                        ilen = iend - istart
                        intron_lens.append(ilen)
                    
                    #onyl take a longest intron!
                    max_intron = max(intron_lens)
                    
                    db_xref = genome[gene]["db_xref"]
                    if isinstance(db_xref, list):
                        db_xref = (",").join(genome[gene]["db_xref"])
                        
                    data = [genome_id, tax_id, taxa["species"], taxa["kingdom"],db_xref, 
                            gene_length, (",").join(list(products)), max_intron]
                    
                    data = list(map(str, data))
                    data = ("\t").join(data)
                    table.write("{}\n".format(data)) 
              
            #taking notes of mammals with no dystrophin gene identified in them
            if not DMD_found:
                data = [genome_id, tax_id, taxa["species"], taxa["kingdom"], "NO DMD Found"]
                data = list(map(str, data))
                data = ("\t").join(data)
                table.write("{}\n".format(data))
                
            
            
            
            
            
            
            
            
            
            
        
            
        
        
        
        
        
        
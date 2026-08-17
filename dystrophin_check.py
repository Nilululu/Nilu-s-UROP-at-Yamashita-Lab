#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 14:23:23 2026

This script will go through all NCBI refeq genomes that we have and filter for mammals. 
Afterwards, it will try to find DMD (dystrophin) gene in each mammals and record gene information
The script takes note of mammals with no DMD gene as well. 
The findings are output in a table named "mammals_dystrophin.txt"

"""

# I need to break this script into multiple functions!

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




def identify_mammal (line):
    """
    takes a line from genomic_directory.csv, gets the metadata of that genome 
    and checks for current mammals 

    Parameters
    ----------
    line : Str : a line from genomic_directory.csv

    Returns
    -------
    None: if it is not a current genome that belons to a mammal
    
    otherwise: 
        loc : Str : location of genomic.gtf
        tax_id : taxonomic id of the genome
        taxa : a dictionary containing taxonomic lineage 
    """
    #parsing a line from genomic_directory
    fields = line.split(",")
    loc = Path(fields[0])
    
    #finding tax_id and status
    meta = get_genome_metadata(loc)
    status = meta[2]
    if status != "current":
        logger.error("not current genome found: {}".format(loc)) #skipping susspended or previous genomes
        return
                  
    tax_id = meta[0]
    taxa = find_taxonomy(tax_id, taxonomy_dict, tax_to_name, {})
    
    #checking for mammals 
    taxa_class = taxa.get("class", "No Class")
    if taxa_class == "mammals": 
        return loc, tax_id, taxa
    
    return 
   
    
   
    
def DMD_parse (gene):
    """
    takes a DMD gene from a genome dictionary created by extract module, 
    and take information of interest out of it

    Parameters
    ----------
    gene : python dict

    Returns
    -------
    info : python list : db_xref, gene_length, max_intron, protein_id
    """
    
    start, end = gene["position"]
    gene_length = end - start
    
    # # do we need this since DMD gene makes dystrophin only? need to check with Romain
    # products = set()
    # for product in gene["products"]:
    #     if product.startswith("dystrophin"):
    #         products.add("dystrophin")
    #     else:
    #         products.add(product)
            
    #finding the max intron on DMD gene
    introns = gene["introns"]
    
    intron_lens = []
    for intron in introns:
        istart, iend = intron
        ilen = iend - istart
        intron_lens.append(ilen)

    max_intron = max(intron_lens)
    #how to find which number intron is the longest one?
    
    #getting gene ids (might have multiple if the gene canbe found in multiple databases)
    db_xref = gene["db_xref"]
    if isinstance(db_xref, list):
        db_xref = (",").join(gene["db_xref"])
    
    #finding the longest transcript:
    max_trc, max_trc_len = None, 0 
    for trc in gene["transcripts"]:
        trc_start, trc_end = gene["transcripts"][trc]["position"]
        trc_len = trc_end - trc_start
        
        if trc_len > max_trc_len:
            max_trc, max_trc_len = trc, trc_len
    
    #can have two protein ids if both automatic and curated annotations were done for the same protein 
    protein_id = gene["transcripts"][max_trc]["protein_id"]
        
    info = [db_xref, gene_length, max_intron, (",").join(list(protein_id))]
    
    return info
        
        
# using genomic dirsctory and makign a table for data I will use
with open ("genomic_directory.csv", 'r') as directory, open("mammals_dystrophin.txt", 'w') as table:
    
    #making table header
    header = "#genome_id tax_id species kingdom db_xref gene_length max_intron protein_id"
    header = ("\t").join(header.split(" ")) 
    table.write("{}\n".format(header))
    
    
    for line in directory:
        
        # will only exist if the genome is a current mammal genome
        mammal = identify_mammal(line)
        if not mammal:
            continue
                
        loc, tax_id, taxa = mammal
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
                
                info = DMD_parse(genome[gene]) 
                    
                #writing data points of interest in the table
                data = [genome_id, tax_id, taxa["species"], taxa["kingdom"]]
                data.extend(info)       
                data = list(map(str, data))
                data = ("\t").join(data)
                table.write("{}\n".format(data)) 
                
            elif genome[gene].get("gene"): # some only have DMD as their gene name and not gene_id attribute 
                if genome[gene].get("gene").upper() == "DMD":
                    
                    if DMD_found: #there shouldn't be two dystrophin gene on the same genome
                        logger.error("two DMD gene found on the same genome, {}".format(genome_id))
                    
                    DMD_found = True
                    
                    info = DMD_parse(genome[gene]) 
                     
                    #writing data points of interest in the table
                    data = [genome_id, tax_id, taxa["species"], taxa["kingdom"]]
                    data.extend(info)       
                    data = list(map(str, data))
                    data = ("\t").join(data)
                    table.write("{}\n".format(data)) 
                
          
        #taking notes of mammals with no dystrophin gene identified in them
        if not DMD_found:
            data = [genome_id, tax_id, taxa["species"], taxa["kingdom"], "NO DMD Found"]
            data = list(map(str, data))
            data = ("\t").join(data)
            logger.error(data)
            
        
        
        
        
        
        
        
        
            
            
        
            
        
        
        
        
        
        
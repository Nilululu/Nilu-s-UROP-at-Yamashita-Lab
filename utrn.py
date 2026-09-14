#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 13:20:05 2026

This script will go through all NCBI refeq genomes that we have and filter for mammals. 
Afterwards, it will try to find UTRN (utrophin) gene in each mammals and record gene information
The script takes note of mammals with no DMD gene as well. 
The findings are output in a table, name TBD

"""

# I need to break this script into multiple functions!

import logging


logging.basicConfig(filename= "log_utrohopin.txt", level = logging.ERROR, force = True)
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
    takes a line from genomic_directory.csv and filters for mammals

    Parameters
    - line : Str : a line from genomic_directory.csv

    Returns
    - None: if it is not a current mammal genome
    - otherwise: 
        loc : Str : location of genomic.gtf
        tax_id : taxonomic id of the genome
        taxa : a dictionary containing taxonomic lineage 
    """
    fields = line.split(",")    #parsing a line from genomic_directory
    loc = Path(fields[0])
    
    meta = get_genome_metadata(loc) #finding tax_id and status
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
   
    
    
def utrn_parse (gene):
    """
    works on a UTRN gene from a genome dictionary created by extract module, 
    and output information of interest

    Parameters
    - gene : python dict

    Returns
    - info : python list : db_xref, gene_length, max_intron, protein_id
    """
    
    start, end = gene["position"]
    gene_length = end - start
    
    introns = gene["introns"]
    
    intron_lens = []
    for intron in introns:
        istart, iend = intron
        ilen = iend - istart
        intron_lens.append(ilen)

    max_intron = max(intron_lens)
    #how to find which number intron is the longest one?
    
    db_xref = gene["db_xref"] #gene ids available in different datasets 
    if isinstance(db_xref, list):
        db_xref = (",").join(gene["db_xref"])
    
    #finding the longest transcript:
    max_trc, max_trc_len = None, 0 
    for trc in gene["transcripts"]:
        trc_start, trc_end = gene["transcripts"][trc]["position"]
        trc_len = trc_end - trc_start
        
        if trc_len > max_trc_len:
            max_trc, max_trc_len = trc, trc_len
     
    protein_id = gene["transcripts"][max_trc]["protein_id"] 
    info = [db_xref, gene_length, max_intron, (",").join(list(protein_id))]
    
    return info
 

       
        
# using genomic dirsctory and making a table for data I will use
with open ("genomic_directory.csv", 'r') as directory, open("mammals_utrn.txt", 'w') as table:
    
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
        
        #making a dict out of gtf file
        genome_id, genome = extract_id_and_genes(loc)
        compute_intron(genome)
        utrn_found = False  #can we have multiple utrn genes?
        
        #searching for utrophin
        for gene in genome:
        
            if gene.upper() == "UTRN":
                if utrn_found: 
                    logger.error("two utrn gene found on the same genome, {}".format(genome_id))
                
                utrn_found = True
                
                info = utrn_parse(genome[gene]) 
                    
                #writing data points of interest in the table
                data = [genome_id, tax_id, taxa["species"], taxa["kingdom"]]
                data.extend(info)       
                data = list(map(str, data))
                data = ("\t").join(data)
                table.write("{}\n".format(data)) 
             
            # I might try to include this block after seeing a similer accurance for utrn
            # elif genome[gene].get("gene"): # some only have DMD as their gene name and not gene_id attribute 
            #     if genome[gene].get("gene").upper() == "UTRN":
                    
            #         if DMD_found: #there shouldn't be two dystrophin gene on the same genome
            #             logger.error("two DMD gene found on the same genome, {}".format(genome_id))
                    
            #         DMD_found = True
                    
            #         info = DMD_parse(genome[gene]) 
                     
            #         #writing data points of interest in the table
            #         data = [genome_id, tax_id, taxa["species"], taxa["kingdom"]]
            #         data.extend(info)       
            #         data = list(map(str, data))
            #         data = ("\t").join(data)
            #         table.write("{}\n".format(data)) 
                
          
        #taking notes of mammals with no dystrophin gene identified in them
        if not utrn_found:
            data = [genome_id, tax_id, taxa["species"], taxa["kingdom"], "NO utrn Found"]
            data = list(map(str, data))
            data = ("\t").join(data)
            logger.error(data)
            
        
        
        
        
        
        
        
        
            
            
        
            
        
        
        
        
        
        
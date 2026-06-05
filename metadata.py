# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from plot_creator import create_hist, create_scatter
import json


def get_GintronInfo(genome, genome_id, threshold, plot = False):
    """
    extracts giants introns and their relative positions with respect
    to the gene they belong to and stores that information in lists

    Parameters
    genome : pyton dict
    genome_id : string
    threshold : int : indicating the cut-off to define a giant intron

    Returns
    -------
    giant_introns : python list : contians lenghts of the giant introns
    in the genome
    distance_from_start : pyhton list : distance between the gene start and the intron start
    distance_from_end : python list : distance between the gene end and the intron end

    """
    
    giant_introns = []
    
    distance_from_start = []  ###
    distance_from_end = []
    
    for gene in genome.values():
        
        gene_start, gene_end= gene["position"]
        
        for item in gene["introns"]:
            if item != set():  #skipping genes with no introns
                start = item[0]
                end = item [1]
                length = end - start
                    
                if length > threshold:
                    giant_introns.append(length)
                
                    start_distance = start - gene_start
                    end_distance = gene_end - end
                    distance_from_start.append(start_distance)
                    distance_from_end.append(end_distance)
            else:
                continue
    if plot:
        if giant_introns:
            create_hist(giant_introns, "giant intron lenghts for {genome_id}")
            create_hist(distance_from_start, "distance from start for {genome_id}")
            create_hist(distance_from_end, "distance from end for {genome_id}")
            create_scatter(distance_from_start, distance_from_end, giant_introns, "", "", "", "z_label")
        else:
            print ("No giant introns for {genome_id} with the current threshold of", threshold)
    return (giant_introns, distance_from_start, distance_from_end)


# I want to make it more efficient by searching for a list of keys

def find_key(key, dictionary):
    if key in dictionary: return dictionary[key]
    for keys, values in dictionary.items():
        if isinstance(values, dict):
            item = find_key(key, values)
            if item is not None:
                return item
            
   

def get_genomeMetadata (gtf_loc):
    """
    parameters: jason_file 
    returns: genome_lenght, num_chromosome, genome_id
    """
    jsonl_file = gtf_loc.parent.parent / "assembly_data_report.jsonl"
    
    with open(jsonl_file, 'r') as report:
        for line in report:
            A = json.loads(line)
            
            taxId = find_key("taxId", A)
            NumChr = find_key("totalNumberOfChromosomes", A)
            Type = find_key("assemblyType", A)
            Status =  find_key("assemblyLevel", A)
            genome_size = find_key("totalSequenceLength", A)
              
    return (taxId, genome_size, Status, NumChr, Type)
  

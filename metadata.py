# -*- coding: utf-8 -*-
"""
a module with functions for extracting genomic metadata in assembly report file as well as 
information about giant introns (needs an already established threshold)
"""

from plot_creator import create_hist, create_scatter
import json


def get_g_intron_info(genome, genome_id, threshold, plot = False):
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

            
def find_key_multiple(keys_of_interest, dictionary, output_dictionary = None):
    """
    finds values of multiple keys inside a nested dictionary 
    and returns a simple dictionary with keys of interest and their values

    Parameters
    ----------
    keys_of_interest : a list
    dictionary : python nested dictionary 
    
    Returns
    -------
    output_dictionary
    """
    
    if output_dictionary is None:
        output_dictionary = {}
        
    for key in keys_of_interest:
        if key in dictionary:
            output_dictionary[key] = dictionary[key]

    for key, value in dictionary.items():
        if isinstance(value, dict):
            output_dictionary = find_key_multiple(keys_of_interest, value, output_dictionary)
    
    return output_dictionary

def get_genome_metadata (gtf_loc):
    """
    parameters: jason_file 
    returns: 
    tax_id, total_sequence_length, assembly_level, assembly_type, 
    numChr, num_scaffolds, num_contigs, scaffold_n50, contig_n50, gc_percent

    """
    jsonl_file = gtf_loc.parent.parent / "assembly_data_report.jsonl"
    with open(jsonl_file, 'r') as report:
        for line in report:
            A = json.loads(line)
            
            keys_of_interest = ["taxId", "totalNumberOfChromosomes", "assemblyStatus", "assemblyType", 
                                "assemblyLevel", "totalSequenceLength", "numberOfContigs", "gcPercent", 
                                "contigN50", "numberOfScaffolds", "scaffoldN50"]
            metadata = find_key_multiple(keys_of_interest, A)
            
            
            tax_num = metadata["taxId"]
            
            try:
                numChr = metadata["totalNumberOfChromosomes"]
            except:
                numChr = "No Data"
            
            assembly_status = metadata["assemblyStatus"]
            assembly_type = metadata["assemblyType"]
            assembly_level =  metadata["assemblyLevel"]
            total_sequence_length = metadata["totalSequenceLength"]
            num_contigs = metadata["numberOfContigs"]
            gc_percent = metadata["gcPercent"]
            contig_n50 = metadata["contigN50"]
            num_scaffolds = metadata["numberOfScaffolds"]
            scaffold_n50 = metadata["scaffoldN50"]
    
    output = [
        tax_num, total_sequence_length, assembly_status, assembly_level, assembly_type, 
        numChr, num_scaffolds, num_contigs, scaffold_n50, contig_n50, gc_percent
        ]
           
    return output


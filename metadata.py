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


# I want to make it more efficient by searching for a list of keys

def find_key(key, dictionary):
    """
    searches in a nested dictioanry until it finds the desired key and returns 
    the values associated with it

    Parameters
    ----------
    key : str
    dictionary : python dict
        a nested dictionary

    Returns
    -------
    item : str
    the values associated with the key
    """
    if key in dictionary: return dictionary[key]
    for keys, values in dictionary.items():
        if isinstance(values, dict):
            item = find_key(key, values)
            if item is not None:
                return item
            
   

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
            
            tax_num = find_key("taxId", A)
            numChr = find_key("totalNumberOfChromosomes", A)
            assembly_type = find_key("assemblyType", A)
            assembly_level =  find_key("assemblyLevel", A)
            total_sequence_length = find_key("totalSequenceLength", A)
            num_contigs = find_key("numberOfContigs", A)
            gc_percent = find_key("gcPercent", A)
            contig_n50 = find_key("contigN50", A)
            num_scaffolds = find_key("numberOfScaffolds", A)
            scaffold_n50 = find_key("scaffoldN50", A)
    
    output = [
        tax_num, total_sequence_length, assembly_level, assembly_type, 
        numChr, num_scaffolds, num_contigs, scaffold_n50, contig_n50, gc_percent
        ]
              
    return output

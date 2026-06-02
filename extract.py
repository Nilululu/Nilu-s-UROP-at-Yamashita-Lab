# -*- coding: utf-8 -*-
"""
This script uses gtf files, extracts the relevent information and 
stores it in a dictionary format
Created on Mon Feb 23 15:41:13 2026 by nilu
"""
from plot_creator import create_hist, create_scatter
import json
import re
import logging
from datetime import datetime

logger = logging.getLogger(__name__)
date = datetime.today().strftime('%Y-%m-%d_%H:%M:%S')
logging.basicConfig(filename="error_parse_log{}.txt".format(date), level=logging.WARNING)

reg_attr = r'(\w+)\s+"([^"]*)"'
def parse_attr_fields(text, reg=reg_attr):
    """ 
    Parse the last field of a gtf ( the attributes fields) and return a dict key: value. if key is repeated return key: [value]

    :param str text: the attribute fields as a string
    :return: the key value dictionnary
    :rtype: dict[str, str]
    """
    results = {}

    for (key, value) in reg.findall(text):
        try:
            this_key = key.strip()
            if this_key in results:
                v = results[this_key]
                if not isinstance(v, list):
                    results[this_key] = [v]
                results[this_key].append(value)
            else:
                results[this_key] = value.replace('"', "").strip() 
        except:
            logging.error("failed to parse line {}".format(text))
            raise
    return results



def extract_genesAndId(gtf_file):
    """
    opens a gtf file and store all the genes, with their respective 
    transcripst, introns, and exons in dict
    
    ----------
    Parameters
    gtf_file : 
        The Gene transfer format (GTF) is a file format used to hold 
        information about gene structure. It is a tab-delimited 
        text format 

    Returns
    genome_id : string
    genes : dictionary
        all the genes and their stored properties that belong the the 
        input genome gtf file

    """

    genes = {}       
    genome_id = None

    with open(gtf_file, "r") as open_gtf:

        for line in open_gtf:

            if line.startswith("#"): # comment or info!

                if line.startswith("#!genome-build-accession"):    
                    # only storing the genome_id from the whole line
                    genome_id = line.strip().split(maxsplit=1)[1] 

                continue
                
            fields = line.strip().split("\t")
            feature = fields[2]
            start = int(fields[3]) - 1 # gtf are 1 bsed to get to 0 based you need to minus 1 here
            end = int(fields[4])
            info_text = fields[8]   

            try:
                attributes = parse_attr_fields(info_text)
            except:
                print(line)
                raise
   
        
            gene_id = attributes["gene_id"]
            #stores the gene_id as a dict with relevent keys
            #if feature == "gene":
            if gene_id not in genes:
                genes[gene_id] = {}
                genes[gene_id]["transcripts"]= {}
                genes[gene_id]["position"]= (start, end)
                genes[gene_id]["introns"]= set()
                genes[gene_id]["strand"]= fields[6]
            
            if feature == "gene":
                genes[gene_id]["position"]= (start, end)
        
            if feature == "exon":
                transcript_id = attributes["transcript_id"]
                # storing exon positions in its appropraite transcript
                #
                if transcript_id not in genes[gene_id]["transcripts"]:
                    genes[gene_id]["transcripts"][transcript_id] = {"exon": [], "intron": []}
                genes[gene_id]["transcripts"][transcript_id]["exon"].append((start, end))
            
   
    return (genome_id, genes)


def compute_intron(genes):
    """
    calculates introns from exon junction  
    
    -----
    parameters: genes : a dictionary 
    returns : genes : a modified dictionary
    
    """
    for gene_id, gene_dico in genes.items():
        
        for transcript_id, transcript_dico in gene_dico["transcripts"].items():
            exons = transcript_dico["exon"]
            if len(exons) <= 1:
                continue #skipping genes with a single exon
            
            else:
                exons = sorted(exons, key=lambda x: x[0])

                for i in range(len(exons)-1):
                    intron_start= int(exons[i][1]) 
                    intron_end= int(exons[i+1][0])
                    transcript_dico["introns"].add((intron_start, intron_end))

                    # depending we may have to modify that like use a set to get unique combination
                    gene_dico["introns"].add((intron_start, intron_end))

 
    return genes
    
   
# plot should be elsewhere
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

def find_key(key, dictionary):
    if key in dictionary: return dictionary[key]
    for keys, values in dictionary.items():
        if isinstance(values, dict):
            item = find_key(key, values)
            if item is not None:
                return item


# this should be elsewehre
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
            
    #       
        
    return (taxId, genome_size, Status, NumChr, Type)
  
def get_biggest_intron(genome):
    """
    extracts the biggest intron within a genome
    Parameters
    genome : pyton dict
    genome_id : string

    Returns
    -------
    max_intron : int
    """
    # here you loose all informations
    # you don't know which gene has the largest intron.

    max_intron = 0

    for gene_id, gene_info in genome.items():
        if intron := (gene_info.get("intron")):
            for (start, end) in intron:
                length = end - start
                    
                if length > max_intron:
                    max_intron = length

    return max_intron
    
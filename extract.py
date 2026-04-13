# -*- coding: utf-8 -*-
"""
This script uses gtf files, extracts the relevent information and 
stores it in a dictionary format
Created on Mon Feb 23 15:41:13 2026

@author: nilof
"""
from plot_creator import create_hist, create_scatter
import typing


def parse_attr_fields(text):
    """ 
    Parse the last field of a gtf ( the attributes fields) and return a dict key: value.

    use the fact that gtf attribute keys are supposed to be in snake case.

    :param str text: the attribute fields as a string
    :return: the key value dictionnary
    :rtype: dict[str, str]
    """
    results = {}
    spt = text.strip().split('";')
    
    #the note part of this field would create problems, so I specified spt with 6 max split
    #gene_id "GUITHDRAFT_162030"; transcript_id "XM_005836996.1"; db_xref "InterPro:IPR007257"; db_xref "JGIDB:Guith1_162030"; gbkey "CDS"; locus_tag "GUITHDRAFT_162030"; note "Subunit of the GINS complex; Psf2; Subunit of the GINS complex.  Psf2"; orig_transcript_id "gnl|WGS:AEIE|GUITHDRAFT_mRNA162030"; product "hypothetical protein"; protein_id "XP_005837053.1"; exon_number "1";
    for element in spt:
        try:
            if element: # test if the string is not empty, may happen
                (key, value) = element.split(maxsplit=1)
                results[key.strip()] = value.replace('"', "").strip() # strip may be useless here not sur
        except:
            print("error")
            print(spt)
            raise
    return results



def extract_genome_info(gtf_file):
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
            start = int(fields[3])
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
                
        
            if feature == "exon":
                transcript_id = attributes["transcript_id"]
                #storing exon positions in its appropraite transcript
                if transcript_id not in genes[gene_id]["transcripts"]:
                    genes[gene_id]["transcripts"][transcript_id] = []
                genes[gene_id]["transcripts"][transcript_id].append((start, end))
            
   
    return (genome_id, genes)


def compute_intron(genes):
    """
    calculates introns from exon junction  
    
    -----
    parameters: genes : a dictionary 
    returns : genes : a modified dictionary
    rtype: dicttionary
    
    """
    for gene in genes.values():
        
        for transcript in gene["transcripts"]:
            if len(gene["transcripts"][transcript]) <= 1:
                continue #skipping genes with a single exon
            
            else:
                exons = sorted(gene["transcripts"][transcript], key=lambda x: x[0])

                for i in range(len(exons)-1):
                    intron_start= int(exons[i][1])
                    intron_end= int(exons[i+1][0])
                    gene["introns"].add((intron_start, intron_end))
 
    return genes
    
   

def get_genome_giant_introns(genome, genome_id, threshold, plot = False):
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

def main(gtf_file, threshold, plot = False):
    """
    combines multiple function in module extract to create a dictionary
    out of the gtf_file relevent information, compute introns for the 
    genome and identify and plot giant introns
    
    PARAMETERS:
    gtf_file 
        The Gene transfer format (GTF) is a file format used to hold 
        information about gene structure. It is a tab-delimited 
        text format 
    threshhold (int)
        indicating the cut-off to define a giant intron
    
    RETURNS:
        genome_id (str)
        genes: python (dict)
        giant Introns (python list)

    """
    (genome_id, genes) = extract_genome_info(gtf_file=gtf_file)
    genes= compute_intron(genes)
    giant_introns_info= get_genome_giant_introns(genes, genome_id, threshold, plot)

    return genome_id, genes, giant_introns_info
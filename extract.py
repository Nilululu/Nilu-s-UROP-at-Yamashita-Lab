# -*- coding: utf-8 -*-
"""
This module has functions for extracting information of interest about genomes and genes 
and their annotated properties from a gtf files and stores it in a dictionary format

Created on Mon Feb 23 15:41:13 2026 by nilu
"""

from datetime import datetime
import re


reg_attr = r'(\w+)\s+"([^"]*)"'
def parse_attr_fields(text, reg=reg_attr):
    """ 
    Parse the last field of a gtf ( the attributes fields) and return a dict key: value. if key is repeated return key: [value]

    :param str text: the attribute fields as a string
    :return: the key value dictionnary
    :rtype: dict[str, str (list)]
    """
    results = {}

    for (key, value) in re.findall(reg, text):
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
            print("failed to parse line {}".format(text))
            raise
    return results



def extract_id_and_genes(gtf_file):
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
        all the genes and their properties of interest

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
            chromosome = fields[0]    
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
            if gene_id not in genes:
                genes[gene_id] = {}
                genes[gene_id]["transcripts"]= {}
                genes[gene_id]["position"]= (start, end)
                genes[gene_id]["introns"]= set()
                genes[gene_id]["products"] = set()
                genes[gene_id]["strand"]= fields[6]
                genes[gene_id]["chr"] = chromosome
                
                
            if feature == "gene":
                genes[gene_id]["position"]= (start, end)
                
                # store the gene product and function and other gene identifier ids there!
                #not all genes have a gene product, which ones have it then?
                
                #NC_058386.1     Gnomon  transcript      26741260        29120869        .       -       .       
                #gene_id "DMD"; transcript_id "XM_045050793.1"; db_xref "GeneID:101084937"; experiment "COORDINATES: polyA evidence [ECO: 0006239]"; gbkey "mRNA"; gene "DMD"; model_evidence "Supporting evidence includes similarity to: 40 Proteins, and 100% coverage of the annotated genomic feature by RNAseq alignments, including 3 samples with support  for all annotated introns"; product "dystrophin, transcript variant X7"; transcript_biotype "mRNA"; 
                # transcripts and exons seems to have this information!
                
            if feature == "exon":
                
                db_xref = attributes.get("db_xref", "No db_xref") #add to gene
                product = attributes.get("product", None)   # should be added to transcript or gene?

                genes[gene_id]["db_xref"] = db_xref
                if product: 
                    genes[gene_id]["products"].add(product) 
                
                transcript_id = attributes["transcript_id"]  
                # storing exon positions in its appropraite transcript

                if transcript_id not in genes[gene_id]["transcripts"]:
                    genes[gene_id]["transcripts"][transcript_id] = {"exons": [], "introns": []}
                    if product:
                        genes[gene_id]["transcripts"][transcript_id]["product"] = product
                        
                        
                genes[gene_id]["transcripts"][transcript_id]["exons"].append((start, end))
            
   
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
            exons = transcript_dico["exons"]
            if len(exons) <= 1:
                continue #skipping genes with a single exon
            
            else:
                exons = sorted(exons, key=lambda x: x[0])

                for i in range(len(exons)-1):
                    intron_start= int(exons[i][1]) 
                    intron_end= int(exons[i+1][0])
                    
                    if intron_start > intron_end:   #to fix the issue of some introns having negetive lengths
                        intron_start, intron_end = intron_end, intron_start
                        
                    transcript_dico["introns"].append((intron_start, intron_end))

                    # depending we may have to modify that like use a set to get unique combination
                    gene_dico["introns"].add((intron_start, intron_end))

 
    return genes
    

    
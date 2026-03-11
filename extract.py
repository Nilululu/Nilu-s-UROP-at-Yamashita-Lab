# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 15:41:13 2026

@author: nilof

This script uses gtf files, extracts the relevent information and 
stores it in the genome and gene classes
"""

import matplotlib.pyplot as plt
import typing

# Style: when defining a variable the "=" sign should have one white space around.
genomes = {}  # a dictionary to stores genomes

def parse_attr_fields(text):
    """ 
    Parse the last field of a gtf ( the attributes fields) and return a dict key: value.

    use the fact that gtf attribute keys are supposed to be in snake case.

    :param str text: the attribute fields as a string
    :return: the key value dictionnary
    :rtype: dict[str, str]
    """
    results = {}
    spt = text.strip().split(";")
    for element in spt:
        if element: # test if the string is not empty may happend.
            (key, value) = element.split(maxsplit=1)
            results[key.strip()] = value.replace('"', "").strip() # strip may be useless here not sure
            # replae is just for the style. 
    return results



def extract_genome_info(gtf_file):
    #openning the file
    # I would recoment you use a context manager while working with file 
    # it may save you some issue in the future!
    # also your variable name are a bit confusing, 
    # I know that is a bit hard when english is not your native languages.
    # but genome_file is okayish it is an approximation gtf_file would be more appropriate.
    # but genome_text is straight up wrong and this is a variable that point to the opened file (called a handle) not the text
    #opening the file 

    # a dict of genes for the genome
    genes = {}       
    genome_id = None

    with open(gtf_file, "r") as open_gtf:

        # I would recommend testing you have the correct line,
        # that way you don't even have to do it in two times.
        for line in open_gtf:
            if line.startswith("#"):# comment or info!
                if line.startswith("#!genome-build-accession"):
                    # also that will extract the genome_id not the full line.
                    # it use split() and strip() to grap the text after the first space whithout the line return.
                    genome_id = line.strip().split(maxsplit=1)[1] 
                continue
                
            fields = line.strip().split("\t")
            feature = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            info_text = fields[8]   
            attributes = parse_attr_fields(info_text)
            #info_text = info_text.strip().strip(";")  # to clena the end of theline ;


        # #skipping other features we don't need 
        #(for some reason no feature passes this block, even exon or gene)
        # if feature != "exon" or "gene":
        #     continue
        
        # you cannot do that at all... you are making you vulnerable to 
        # the order may change between file and event inside a file.
        # that is why I was suggesitng using a funciton that parse that.
        #info text contains many parameters we don't need
        #we only store gene and transcripts id for each gene and exon
        #info_list= info_text.split("; ")
        #gene_id = info_list[0].split(" ")[1]
        #transcript_id = info_list[1].split(" ")[1]
        
        gene_id = attributes["gene_id"]

        #stores the gene_id as a dict with relevent keys
        if feature == "gene":
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
        
    # so the nice thing with the way I showing you to read file is that you don't have to remeber to close them.
    # they close automatically as soon as you end the matching identation.
    return (genome_id, genes)
# I don't think that belongs to this function 
# you should make a new one.
def compute_intron(genes):
    #calculating introns from exon junction    
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
    # here you only reading one file 
    # if you want to keep track of multiple genome using a dictionany is a good idea, 
    # but it should not be here
    #genomes[genome_id]= genes
    
    #there is a problem here that should be fixed later
    #does this function creates a new genomes dict everytime or reuses the same one
    #return genomes

# by splitting it like I do it may seems unecessary but it will make the code easier to follow and write when we will have
# to write the parser for everything.

# so you may want to make a main function now()
# this will be very convenient if we want to use that in a CLI or to call it from outside (as a module)

def main(gtf_file, threshold):
    (genome_id, genes) = extract_genome_info(gtf_file=gtf_file)
    get_genome_giant_introns(genes, genome_id, threshold)

# try to be consistant with your naming convention
# is it genes or genome?
def get_genome_giant_introns(genome, genome_id, threshold):
    
    giant_introns = []
    
    distance_from_start = []  ###
    distance_from_end = []
    
    for gene in genome.values():
        
        gene_start, gene_end= gene["position"]  ###
        
        for item in gene["introns"]:
            if item != set():
                start = item[0]
                end = item [1]
                length = end - start
                    
                if length > threshold:
                    giant_introns.append(length)
                
                    ###
                    relative_start = start - gene_start
                    end_distance = gene_end - end
                    distance_from_start.append(relative_start)
                    distance_from_end.append(end_distance)
            else:
                continue
    # ok let's talk about layout next time!
    fig, ax = plt.subplots(1)  #for intron lenghts  
    fig1, ax1 = plt.subplots(1)  #for introns distance from gene start
    fig2, ax2 = plt.subplots(1)  #for introns distance from gene end
    fig3, ax3 = plt.subplots(1)  #for all the above
    
    ax.hist(giant_introns, int(max(giant_introns)/1000))
    ax.set_title("Giant Introns Lenghts")
    
    ax1.hist(distance_from_start, 100)
    ax1.set_title("Giant Introns start relative to their gene")
    
    ax2.hist(distance_from_end, 100) 
    ax2.set_title("Giant Intron ends distance from gene ends")
    
    sc= ax3.scatter(distance_from_start, distance_from_end, c=giant_introns, cmap="viridis")
    ax3.set_xlabel("Introns start relative to their gene start")
    ax3.set_ylabel("Intron ends distance from gene ends")

    cbar = plt.colorbar(sc)
    cbar.set_label('Intron Lengths')

    plt.title(genome_id)
   
    plt.tight_layout()
    plt.show()

    
   
    return giant_introns



    
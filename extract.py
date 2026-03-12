# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 15:41:13 2026

@author: nilof

This script uses gtf files, extracts the relevent information and 
stores it in the genome and gene classes
"""

import matplotlib.pyplot as plt


genomes= {}  # a dictionary to stores genomes

def extract_genome_info(genome_file):
    #openning the file
    genome_text = open(genome_file, "r")  #reading the file
    genome_id=""
    
    
    #get the genome ID from line 3
    for i, line in enumerate(genome_text, start= 1):
        if i==3:
            genome_id= line
            break
    
    # a dict of genes for the genome
    genes = {}   
    
    #extracting genes and exons
    for line in genome_text:
        
        #skipping informational lines
        if line.startswith("#"):
            continue
        
        fields = line.strip().split("\t")
        feature = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        info_text = fields[8]   
        info_text = info_text.strip().strip(";")  # to clena the end of theline ;
        
        

        # #skipping other features we don't need 
        #(for some reason no feature passes this block, even exon or gene)
        # if feature != "exon" or "gene":
        #     continue
        
        #info text contains many parameters we don't need
        #we only store gene and transcripts id for each gene and exon
        info_list= info_text.split("; ")
        gene_id = info_list[0].split(" ")[1]
        transcript_id = info_list[1].split(" ")[1]
        
        
        #stores the gene_id as a dict with relevent keys
        if feature == "gene":
            if gene_id not in genes:
                genes[gene_id] = {}
                genes[gene_id]["transcripts"]= {}
                genes[gene_id]["position"]= (start, end)
                genes[gene_id]["introns"]= set()
                
            
            
        
        if feature == "exon":
            
            #storing exon positions in its appropraite transcript
            if transcript_id not in genes[gene_id]["transcripts"]:
                genes[gene_id]["transcripts"][transcript_id] = []
            
            genes[gene_id]["transcripts"][transcript_id].append((start, end))
        
           
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
   
    genomes[genome_id]= genes
    
    #there is a problem here that should be fixed later
    #does this function creates a new genomes dict everytime or reuses the same one
    return genomes




def get_genome_giant_introns(genome, threshold):
    
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
   
    plt.tight_layout()
    plt.show()

    
    # need to work on what to return
    return giant_introns



    
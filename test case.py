# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 15:41:16 2026 by nilu

tests the modulaes "extract" and "plot_creator" 
uses Guilt genome: the smallest genome among eaukaryotes 
"""

from extract import extract_genesAndId, compute_intron, get_GintronInfo, get_genomeMetadata
from pathlib import Path
from plot_creator import create_hist, create_scatter, create_2d_scatter, get_style
import numpy as np
import taxonomy 
import matplotlib.pyplot as plt
import math
#### incorporate the taxonomy module 


#for using the script in cluster
nodes_file = "nodes.dmp"
names_file = "names.dmp"

taxonomy_dict = taxonomy.generate_taxonomy_dict(nodes_file)
tax_to_name = taxonomy.generate_tax_to_name(names_file)

genomes = dict()
Gintron_treshold = 1000

allGenomes_Gintron = []
allGenomes_GintronStart = []
allGenomes_GintronEnd = []

GenomesGintrons_dict = {}
genomeId_to_taxonomy = {}
genome_Metadata = {}

#  file containing all gtf locations
#!!!!! 
csv_file = "all_gtf_downloaded_find.txt"  

# csv_file = "genomic_directory.csv"    #used on local drive

with open(csv_file, mode = 'r') as directory:
    lines = directory.readlines()

    for line in lines[:200]:
       
        try:
            
            # name, gtf_loc = line.split(",")
            gtf_loc = line.strip()    #used on cluster  #!!!!!!
            gtf_loc = Path (gtf_loc.strip())
            
            genomeId, genes = extract_genesAndId(gtf_loc)
            compute_intron(genes)
            
            Gintron, GintronStart, GintronEnd = get_GintronInfo(genes, genomeId, Gintron_treshold)
            
            # print(taxId)
            taxNum, genome_size, status, numChr, genome_type= get_genomeMetadata(gtf_loc)
            taxId = taxonomy.find_taxonomy(taxNum, taxonomy_dict, tax_to_name)
            
            genomes[genomeId]= genes
            
            
            if Gintron:
                allGenomes_Gintron.extend(Gintron)
                allGenomes_GintronStart.extend(GintronStart)
                allGenomes_GintronEnd.extend(GintronEnd)
            
            GenomesGintrons_dict[genomeId] = {}
            GenomesGintrons_dict[genomeId]["Gintrons"] = np.column_stack((Gintron, GintronStart, GintronEnd))
            
            genome_Metadata[genomeId] = {}
            genome_Metadata[genomeId]["size"] = int(genome_size)
            genome_Metadata[genomeId]["Num_Chr"] = int(numChr)
            genome_Metadata[genomeId]["taxId"]  = taxId
            genome_Metadata[genomeId]["genome_type"] = genome_type
       
        except Exception as e:
            print("failed to extract the information for genome:", line)
            print("due to this error:", e)
            raise 
        
        
allGenomes_GintronMatrix = np.column_stack((allGenomes_Gintron, allGenomes_GintronStart, allGenomes_GintronEnd))


fig, ax = plt.subplots(1)
fig1, ax1 = plt.subplots(1)

_mapping = dict()
_kingdom_count = dict()
_type_count = dict()

for genome in GenomesGintrons_dict:
    # print(genome)
    size = genome_Metadata[genome]["size"]
    
    max_intron = 0
    
    # print ("size", GenomesGintrons_dict[genome]["Gintrons"].size)
    
    
    if GenomesGintrons_dict[genome]["Gintrons"].size > 0:
        Gintron_f = GenomesGintrons_dict[genome]["Gintrons"]
        max_intron = Gintron_f[:,0].max()
    
    else:
        
        try: 
            for gene in genomes[genome]:
                # print(genomes[genome][gene])
                for item in genomes[genome][gene]['introns']:
                    # print("item", item)
                    if item != set():  #skipping genes with no introns
                        start = item[0]
                        end = item [1]
                        length = end - start
                        
                        # print("length", length)
                            
                        if length > max_intron:
                            max_intron = length
                        
                    else:
                        continue
            # print("success", genome)
        
        except Exception as e:
            print(item)
            print("failed attempt to extract introns", genome, e)
            print(genomes[genome][gene].keys())
            # raise

        
        
   
    kingdom = genome_Metadata[genome]["taxId"]["kingdom"]
    genome_type = genome_Metadata[genome]["genome_type"]
    
    try:
        _mapping = get_style(kingdom)
    
    except Exception as e:
        print("failed to get style for:", kingdom, genome, e)
        print(_mapping)
        
    style = _mapping[kingdom]
    if kingdom not in _kingdom_count:
        _kingdom_count[kingdom] = 1
    
    else:
        _kingdom_count[kingdom] += 1
        
    if genome_type not in _type_count:
        _type_count[genome_type] = 1
    
    else:
        _type_count[genome_type] +=1
        
    size = genome_Metadata[genome]["size"]
    
    # print(size)
    # print(size, max_intron, style['color'], style['marker'] )
    ax.scatter((size*0.001), (max_intron*0.001), color = style['color'], marker = style['marker'])
    ax.set_title(f"genome size vs max intron (scale = kbp), labels: {_mapping}")
    ax1.scatter(math.log10(size*0.001), math.log10(max_intron*0.001), color = style['color'], marker = style['marker'] )
    ax1.set_title("genome size vs max intron (scale = log bp)")

plt.show()  
print(_mapping)   
    
fig2, ax2 = plt.subplots(1)
print(_kingdom_count)
ax2.pie(_kingdom_count.values(), labels = _kingdom_count.keys())

fig3, ax3 = plt.subplots(1)
ax3.pie(_type_count.values(), labels = _type_count.keys())
    
ax.savefig("genome_size_vs_max_intron")
ax1.savefig("genome_size_vs_max_intron_log")
ax2.savefig("kingdom_distribution")
ax3.savefig("type_distribution") 
    
    
    

# #### Filteration for all questions 
# filter_name = "Mycota"
# filter_rank = "kingdom"
# filteration = False

# ### Question 1 = What is the distribution of number of Gintrons in each genome




    
# if filteration == True:
    
#     filt_numGintrons = []
#     filter_Gintrons = []
#     filter_GintronStart = []
#     filter_GintronEnd = []
#     filt_GenomeLength = []
#     filt_Genomes_numChr = []
#     filt_max_Gintrons = []
#     filt_mean_Gintrons = []
#     filt_quantile_list = []
#     quantile = 0.99
    
#     for item in GenomesGintrons_dict:
#         if GenomesGintrons_dict[item]["taxId"][filter_rank] == filter_name:
    
#             numGintrons = GenomesGintrons_dict[item]["Gintrons"].shape[0]
#             filt_numGintrons.append(numGintrons)
            
#             Gintron_f = GenomesGintrons_dict[item]["Gintrons"]
#             filter_Gintrons.extend(Gintron_f[:,0])
#             filter_GintronStart.extend(Gintron_f[:,1])
#             filter_GintronEnd.extend(Gintron_f[:,2])
#             filt_GenomeLength.append(GenomesGintrons_dict [item]["GenomeLength"])
#             filt_Genomes_numChr.append(GenomesGintrons_dict[item]["Num_Chr"])
            
#             try: 
#                 max_Gintron = Gintron_f[:,0].max()
#                 mean_Gintron = Gintron_f[:,0].mean()
#                 quantile_Gintron = np.quantile(Gintron_f[:,0], quantile)
#             except:
#                 max_Gintron = 0
#                 mean_Gintron = 0
#                 quantile_Gintron = 0 
            
            
#             filt_max_Gintrons.append(max_Gintron)
#             filt_mean_Gintrons.append(mean_Gintron)
#             filt_quantile_list.append(quantile_Gintron)
            
#     if filt_numGintrons:
#         create_hist(filt_numGintrons, f"Q1: distribution of number of giant_introns in genomes (filter : {filter_name} {filter_rank}, Gintron threshold: {Gintron_treshold}kb)")
        
#         ### Question 2 = for all Gintrons, what is their length distribution 
        
#         create_hist(filter_Gintrons, f"Q2: giant introns lenght distribution (filter : {filter_name} {filter_rank},Gintron threshold: {Gintron_treshold}kb)")
          
        
#         ### Question 3 = Is there a corelation between Gintron length and its distance from two ends of the gene
        
#         Q3_title  = f"Q3: Gintrons lenght vs Distance from both end of a gene (filter : {filter_name} {filter_rank},Gintron threshold: {Gintron_treshold}kb)"
#         Q3_xlabel = "Distance from start"
#         Q3_ylabel = "Distance from end"
#         Q3_zlabel = "Gintron Lenght"
        
#         print(filter_Gintrons, filter_GintronStart,filter_GintronEnd)
#         create_scatter(filter_GintronStart, filter_GintronEnd, filter_Gintrons, Q3_title, Q3_xlabel, Q3_ylabel, Q3_zlabel)
        
        
        
#         ## Question 4 = Is there a correlation between genome lenght&number of chromosomes and number of Gintrons, max Gintron, mean GIntron, different quantiles of Gintrons
        
#         create_scatter(filt_GenomeLength, filt_numGintrons, filt_Genomes_numChr, f"Q4.a genome len vs num Gintrons vs num chr (filter : {filter_name} {filter_rank}, Gintron threshold: {Gintron_treshold}kb)", "genome len", "num Gintrons" , "num chr")
#         create_scatter(filt_GenomeLength , filt_max_Gintrons, filt_Genomes_numChr,  f"Q4b.b Genome len vs biggest Gintron vs num chr (filter : {filter_name} {filter_rank}, Gintron threshold: {Gintron_treshold}kb)", "genome len", "max Gintron" , "num chr")
#         create_scatter(filt_GenomeLength , filt_mean_Gintrons, filt_Genomes_numChr,  f"Q4.c Genome len vs mean Gintron vs num chr wiht (filter : {filter_name} {filter_rank}, Gintron threshold: {Gintron_treshold}kb)", "genome len", "max Gintron" , "num chr")
#         create_scatter(filt_GenomeLength , filt_quantile_list, filt_Genomes_numChr,  f"Q4.d Genome len vs {quantile}-th quantile Gintron vs num chr (filter : {filter_name} {filter_rank}, Gintron threshold: {Gintron_treshold}kb)", "genome len", f"{quantile}-th quantile Gintron" , "num chr")
        
        
#     else:
#         print("No genomes found with the current filter") 
# if filteration == False:
#     allGenomes_numGintrons = []
    
#     for item in GenomesGintrons_dict:
        
#         numGintrons = GenomesGintrons_dict[item]["Gintrons"].shape[0]
#         allGenomes_numGintrons.append(numGintrons)
#     create_hist(allGenomes_numGintrons, f"Q1: distribution of number of giant_introns in genomes (Gintron threshold: {Gintron_treshold} bp)")


#     ### Question 2 = for all Gintrons, what is their length distribution 
    
#     create_hist(allGenomes_Gintron, f"Q2: giant introns lenght distribution (Gintron threshold: {Gintron_treshold}bp)")

    
    
    
    
#     ### Question 3 = Is there a corelation between Gintron length and its distance from two ends of the gene
    
#     Q3_title  = f"Q3: Gintrons lenght vs Distance from both end of a gene (Gintron threshold: {Gintron_treshold}bp)"
#     Q3_xlabel = "Distance from start"
#     Q3_ylabel = "Distance from end"
#     Q3_zlabel = "Gintron Lenght"
    
#     create_scatter(allGenomes_GintronMatrix[:,1], allGenomes_GintronMatrix[:,2], allGenomes_GintronMatrix[:,0], Q3_title, Q3_xlabel, Q3_ylabel, Q3_zlabel)
    
    
    
#     ### Question 4 = Is there a correlation between genome lenght&number of chromosomes and number of Gintrons, max Gintron, mean GIntron, different quantiles of Gintrons
    
#     quantile = 80
    
#     allGenomes_numGintrons
#     allGenomes_GenomeLength = []
#     allGenomes_numChr = []
#     max_Gintrons = []
#     mean_Gintrons = []
#     quantile_list = []
    
#     for key in GenomesGintrons_dict:
        
#         GenomeLength = GenomesGintrons_dict [key]["GenomeLength"]
#         numChr = GenomesGintrons_dict[key]["Num_Chr"]
       
#         try: 
#             max_Gintron = GenomesGintrons_dict[key]["Gintrons"].max()
#             mean_Gintron = GenomesGintrons_dict[key]["Gintrons"].mean()
#             quantile_Gintron = np.quantile(np.array(GenomesGintrons_dict[key]["Gintrons"]), 0.75)
#         except:
#             max_Gintron = 0
#             mean_Gintron = 0
#             quantile_Gintron = 0 
            
#         allGenomes_GenomeLength.append(GenomeLength)
#         allGenomes_numChr.append(numChr)
#         max_Gintrons.append(max_Gintron)
#         mean_Gintrons.append(mean_Gintron)
#         quantile_list.append(quantile_Gintron)
    
#     create_scatter(allGenomes_GenomeLength, allGenomes_numGintrons, allGenomes_numChr, f"Q4.a genome len vs num Gintrons vs num chr (Gintron threshold: {Gintron_treshold}bp)", "genome len", "num Gintrons" , "num chr")
#     create_scatter(allGenomes_GenomeLength , max_Gintrons, allGenomes_numChr,  f"Q4b.b Genome len vs biggest Gintron vs num chr (Gintron threshold: {Gintron_treshold}bp)", "genome len", "max Gintron" , "num chr")
#     create_scatter(allGenomes_GenomeLength , mean_Gintrons, allGenomes_numChr,  f"Q4.c Genome len vs mean Gintron vs num chr wiht (Gintron threshold: {Gintron_treshold}bp)", "genome len", "max Gintron" , "num chr")
#     create_scatter(allGenomes_GenomeLength , quantile_list, allGenomes_numChr,  f"Q4.d Genome len vs {quantile}-th quantile Gintron vs num chr (Gintron threshold: {Gintron_treshold}bp)", "genome len", f"{quantile}-th quantile Gintron" , "num chr")
    
  
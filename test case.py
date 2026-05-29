# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 15:41:16 2026 by nilu

tests the modulaes "extract" and "plot_creator" 
uses Guilt genome: the smallest genome among eaukaryotes 
"""

from extract import extract_genesAndId, compute_intron, get_GintronInfo, get_taxId, get_genomeMetadata, get_biggest_intron
from pathlib import Path
from plot_creator import create_hist, create_scatter, create_2d_scatter, get_style
import numpy as np
import taxonomy 
import matplotlib.pyplot as plt

#### incorporate the taxonomy module 

# for when I try the script in my local drive
# nodes_file = Path(r"C:\Users\nilof\OneDrive\Documents\Python Scripts\new_taxdump\nodes.dmp")
# names_file = Path(r"C:\Users\nilof\OneDrive\Documents\Python Scripts\new_taxdump\names.dmp")

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

#  file containing all gtf locations
csv_file = "all_gtf_downloaded_find.txt"

with open(csv_file, mode = 'r') as directory:
    lines = directory.readlines()

    for line in lines:
       
        name, gtf_loc = line.split(",")
        gtf_loc = Path (gtf_loc.strip())
        
        genomeId, genes = extract_genesAndId(gtf_loc)
        compute_intron(genes)
        
        Gintron, GintronStart, GintronEnd = get_GintronInfo(genes, genomeId, Gintron_treshold)
        taxNum = get_taxId(gtf_loc)
        taxId = taxonomy.find_taxonomy(taxNum, taxonomy_dict, tax_to_name)
        
        # print(taxId)
        genomeMetadata = get_genomeMetadata(gtf_loc)

        genomeId_to_taxonomy = taxId
        genomes[genomeId]= genes
        
        
        if Gintron:
            allGenomes_Gintron.extend(Gintron)
            allGenomes_GintronStart.extend(GintronStart)
            allGenomes_GintronEnd.extend(GintronEnd)
        
        GenomesGintrons_dict[genomeId] = {}
        GenomesGintrons_dict[genomeId]["Gintrons"] = np.column_stack((Gintron, GintronStart, GintronEnd))
        GenomesGintrons_dict[genomeId]["GenomeLength"] = genomeMetadata[0]
        GenomesGintrons_dict[genomeId]["Num_Chr"] = genomeMetadata[1]
        GenomesGintrons_dict[genomeId]["taxId"]  = taxId

allGenomes_GintronMatrix = np.column_stack((allGenomes_Gintron, allGenomes_GintronStart, allGenomes_GintronEnd))


fig, ax = plt.subplots(1)
# fig1, ax1 = plt.subplots(1)

_mapping = dict()

for genome in GenomesGintrons_dict:
    print(genome)
    size = GenomesGintrons_dict[genome]["GenomeLength"]
    
    max_intron = 0
    
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

        
        
   
    kingdom = GenomesGintrons_dict[genome]["taxId"]["kingdom"]
    _mapping = get_style(kingdom)
    
    style = _mapping[kingdom]
    
    # print(size, max_intron, style['color'], style['marker'] )
    ax.scatter(size, max_intron, color = style['color'], marker = style['marker'])
    
print(_mapping)   
    

    
    
    
    
    
    
    

# #### Filteration for all questions 
# filter_name = "Mycota"
# filter_rank = "kingdom"
# filteration = True

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
        
#         # Q5 what does the above questions look like for different filterations when you take kingdoms and their sub catagpries in mind?
#         # I need a filteration function !!!
        
#     else:
#         print("No genomes found with the current filter") 
# if filteration == False:
#     allGenomes_numGintrons = []
    
#     for item in GenomesGintrons_dict:
        
#         numGintrons = GenomesGintrons_dict[item]["Gintrons"].shape[0]
#         allGenomes_numGintrons.append(numGintrons)
#     create_hist(allGenomes_numGintrons, f"Q1: distribution of number of giant_introns in genomes (Gintron threshold: {Gintron_treshold}kb)")


#     ### Question 2 = for all Gintrons, what is their length distribution 
    
#     create_hist(allGenomes_Gintron, f"Q2: giant introns lenght distribution (Gintron threshold: {Gintron_treshold}kb)")

    
    
    
    
#     ### Question 3 = Is there a corelation between Gintron length and its distance from two ends of the gene
    
#     Q3_title  = f"Q3: Gintrons lenght vs Distance from both end of a gene (Gintron threshold: {Gintron_treshold}kb)"
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
    
#     create_scatter(allGenomes_GenomeLength, allGenomes_numGintrons, allGenomes_numChr, f"Q4.a genome len vs num Gintrons vs num chr (Gintron threshold: {Gintron_treshold}kb)", "genome len", "num Gintrons" , "num chr")
#     create_scatter(allGenomes_GenomeLength , max_Gintrons, allGenomes_numChr,  f"Q4b.b Genome len vs biggest Gintron vs num chr (Gintron threshold: {Gintron_treshold}kb)", "genome len", "max Gintron" , "num chr")
#     create_scatter(allGenomes_GenomeLength , mean_Gintrons, allGenomes_numChr,  f"Q4.c Genome len vs mean Gintron vs num chr wiht (Gintron threshold: {Gintron_treshold}kb)", "genome len", "max Gintron" , "num chr")
#     create_scatter(allGenomes_GenomeLength , quantile_list, allGenomes_numChr,  f"Q4.d Genome len vs {quantile}-th quantile Gintron vs num chr (Gintron threshold: {Gintron_treshold}kb)", "genome len", f"{quantile}-th quantile Gintron" , "num chr")
    
    # Q5 what does the above questions look like for different filterations when you take kingdoms and their sub catagpries in mind?
    # I need a filteration function !!! 
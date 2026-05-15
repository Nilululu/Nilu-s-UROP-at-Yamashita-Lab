# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 15:41:16 2026 by nilu

tests the modulaes "extract" and "plot_creator" 
uses Guilt genome: the smallest genome among eaukaryotes 
"""

from extract import extract_genesAndId, compute_intron, get_GintronInfo, get_taxId, get_genomeMetadata
from pathlib import Path
from plot_creator import create_hist, create_scatter, create_2d_scatter
import numpy as np




genomes = dict()
Gintron_treshold = 1000

allGenomes_Gintron = []
allGenomes_GintronStart = []
allGenomes_GintronEnd = []

GenomesGintrons_dict = {}


#  file containing all gtf locations
csv_file = "genomic_directory.csv"

with open(csv_file, mode = 'r') as directory:
    lines = directory.readlines()

    for line in lines:
       
        name, gtf_loc = line.split(",")
        gtf_loc = Path (gtf_loc.strip())
        
        genomeId, genes = extract_genesAndId(gtf_loc)
        compute_intron(genes)
        
        Gintron, GintronStart, GintronEnd = get_GintronInfo(genes, genomeId, Gintron_treshold)
        taxId = get_taxId(gtf_loc)
        # genomeId, genes, giant_introns_info, taxId = main(Path(gtf_loc.strip()), Gintron_treshold)
        genomes[genomeId]= genes
        # giant_introns_len, start_distance, end_distance = giant_introns_info
        
        # report_path = gtf_loc.replace("genomic.gtf", "sequence_report.jsonl").strip()
        genomeMetadata= get_genomeMetadata(gtf_loc)
        
        if Gintron:
            allGenomes_Gintron.extend(Gintron)
            allGenomes_GintronStart.extend(GintronStart)
            allGenomes_GintronEnd.extend(GintronEnd)
            
        GenomesGintrons_dict[genomeId] = {}
        GenomesGintrons_dict[genomeId]["Gintrons"] = np.column_stack((Gintron, GintronStart, GintronEnd))
        GenomesGintrons_dict[genomeId]["GenomeLength"] = genomeMetadata[0]
        GenomesGintrons_dict[genomeId]["Num_Chr"] = genomeMetadata[1]

allGenomes_GintronMatrix = np.column_stack((allGenomes_Gintron, allGenomes_GintronStart, allGenomes_GintronEnd))



### Question 1 = What is the distribution of number of Gintrons in each genome

allGenomes_numGintrons = []

for item in GenomesGintrons_dict:
    
    numGintrons = GenomesGintrons_dict[item]["Gintrons"].shape[0]
    allGenomes_numGintrons.append(numGintrons)

create_hist(allGenomes_numGintrons, f"distribution of number of giant_introns in genomes wiht {Gintron_treshold}kb treshold")




### Question 2 = for all Gintrons, what is their length distribution 

create_hist(allGenomes_Gintron, "giant introns lenght distribution wiht {Gintron_treshold}kb as treshold")





### Question 3 = Is there a corelation between Gintron length and its distance from two ends of the gene

Q3_title  = f"Gintrons lenght vs Distance from both end of a gene, {Gintron_treshold}kb as treshold"
Q3_xlabel = "Distance from start"
Q3_ylabel = "Distance from end"
Q3_zlabel = "Gintron Lenght"

create_scatter(allGenomes_GintronMatrix[:,1], allGenomes_GintronMatrix[:,2], allGenomes_GintronMatrix[:,0], Q3_title, Q3_xlabel, Q3_ylabel, Q3_zlabel)



### Question 4 = Is there a correlation between genome lenght&number of chromosomes and number of Gintrons, max Gintron, mean GIntron, different quantiles of Gintrons

quantile = 75

allGenomes_numGintrons
allGenomes_GenomeLength = []
allGenomes_numChr = []
max_Gintrons = []
mean_Gintrons = []
quantile_list = []

for key in GenomesGintrons_dict:
    
    GenomeLength = GenomesGintrons_dict [key]["GenomeLength"]
    numChr = GenomesGintrons_dict[key]["Num_Chr"]
   
    try: 
        max_Gintron = GenomesGintrons_dict[key]["Gintrons"].max()
        mean_Gintron = GenomesGintrons_dict[key]["Gintrons"].mean()
        quantile_Gintron = np.quantile(np.array(GenomesGintrons_dict[key]["Gintrons"]), 0.75)
    except:
        max_Gintron = 0
        mean_Gintron = 0
        quantile_Gintron = 0 
        
    allGenomes_GenomeLength.append(GenomeLength)
    allGenomes_numChr.append(numChr)
    max_Gintrons.append(max_Gintron)
    mean_Gintrons.append(mean_Gintron)
    quantile_list.append(quantile_Gintron)

create_scatter(allGenomes_GenomeLength, allGenomes_numGintrons, allGenomes_numChr, f"genome len vs num Gintrons vs num chr wiht {Gintron_treshold}kb as treshold", "genome len", "num Gintrons" , "num chr")
create_scatter(allGenomes_GenomeLength , max_Gintrons, allGenomes_numChr,  f"Genome len vs biggest Gintron vs num chr wiht {Gintron_treshold}kb as treshold", "genome len", "max Gintron" , "num chr")
create_scatter(allGenomes_GenomeLength , mean_Gintrons, allGenomes_numChr,  f"Genome len vs mean Gintron vs num chr wiht {Gintron_treshold}kb as treshold", "genome len", "max Gintron" , "num chr")
create_scatter(allGenomes_GenomeLength , mean_Gintrons, allGenomes_numChr,  f"Genome len vs mean Gintron vs num chr wiht {Gintron_treshold}kb as treshold", "genome len", "mean Gintron" , "num chr")
create_scatter(allGenomes_GenomeLength , quantile_list, allGenomes_numChr,  f"Genome len vs {quantile}-th quantile Gintron vs num chr wiht {Gintron_treshold}kb as treshold", "genome len", f"{quantile}-th quantile Gintron" , "num chr")
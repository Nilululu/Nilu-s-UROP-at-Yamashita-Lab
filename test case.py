# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 15:41:16 2026

@author: nilof
"""
"""
tests the modulaes "extract" and "plot_creator" 
uses Guilt genome: the smallest genome among eaukaryotes 
"""

from extract import extract_genome_info
from extract import get_genome_giant_introns
from extract import main
from extract import compute_intron
from pathlib import Path
from plot_creator import create_hist, create_scatter
import matplotlib.pyplot as plt
import numpy as np

# for each line open the location given for each genome's gtf file
#extract information from gtf
#add the extracted information to a dictionary
#use that dictionary to ask big picture questions 
genomes = dict()
giant_introns_list = []
start_distance_list = []
end_distance_list = []
giant_introns_dict = {}

# read the cvs 
csv_file = "genomic_directory.csv"

with open(csv_file, mode = 'r') as directory:
    lines = directory.readlines()

    for line in lines:
       
        name, gtf_loc = line.split(",")
        genome_id, genome, giant_introns_info = main(Path(gtf_loc.strip()), 1000)
        genomes[genome_id]= genome
        giant_introns, start_distance, end_distance = giant_introns_info
        giant_introns_list.append(giant_introns)
        start_distance_list.append(start_distance)
        end_distance_list.append(end_distance)
        giant_introns_dict[genome_id] = giant_introns


#number of giant introns inside each gene 
num_g_introns_list = []
for items in giant_introns_dict:
    
    num_g_introns = len(giant_introns_dict[items])
    num_g_introns_list.append(num_g_introns)

print(num_g_introns_list)


#asking number of giant introns/num introns for genomes
intron_ratio= []
introns = []

for id_ in giant_introns_dict:
    num_g_introns = len(giant_introns_dict[id_])

for item in genomes.values():
    
    genome_introns = 0
    for gene in item.values():
        num_introns = len(gene['introns'])
        genome_introns = genome_introns + num_introns
    introns.append(genome_introns)


a = np.array(num_g_introns_list)
b = np.array(introns)

intron_ratio = a/b
intron_ratio = list(intron_ratio)

print("intron ratio", intron_ratio)
print("mean", np.mean(intron_ratio)*100)
print("std", np.std(intron_ratio)*100)

create_hist(giant_introns_list, "giant introns lenght distribution")
create_hist(introns,"number of introns per genome distribution") 
create_hist(intron_ratio, "ration of giant introns within all introns per genome distribution")

print("dimentions", len(giant_introns_list), len(start_distance_list), len(end_distance_list))

print(start_distance_list, end_distance_list, giant_introns_list)
# create_scatter(start_distance_list, end_distance_list, giant_introns_list, "title", "start distance", "end distance", "giant intron lenght")

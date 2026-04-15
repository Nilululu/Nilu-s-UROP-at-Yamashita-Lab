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

treshold = 1000

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
        genome_id, genome, giant_introns_info = main(Path(gtf_loc.strip()), treshold)
        genomes[genome_id]= genome
        giant_introns_len, start_distance, end_distance = giant_introns_info
        if giant_introns_len:
            giant_introns_list.extend(giant_introns_len)
            start_distance_list.extend(start_distance)
            end_distance_list.extend(end_distance)
        giant_introns_dict[genome_id] = np.column_stack((giant_introns_len, start_distance, end_distance))
        


giant_introns = np.column_stack((giant_introns_list, start_distance_list, end_distance_list))

#number of giant introns inside each genome 
num_g_introns = []

for items in giant_introns_dict:
    
    g_introns = giant_introns_dict[items].shape[0]
    
    num_g_introns.append(g_introns)

create_hist(num_g_introns, f"distribution of number of giant_introns in genomes wiht {treshold}kb treshold")



create_hist(giant_introns_list, "giant introns lenght distribution")


create_scatter(giant_introns[:,0], giant_introns[:,1], giant_introns[:,2], "title", "start distance", "end distance", "giant intron lenght")

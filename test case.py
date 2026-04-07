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
import pathlib


# for each line open the location given for each genome's gtf file
#extract information from gtf
#add the extracted information to a dictionary
#use that dictionary to ask big picture questions 
genomes ={}

# read the cvs 
csv_file = "genomic_directory.csv"

with open(csv_file, mode = 'r') as directory:
    lines = directory.readlines()

    for line in lines:
       
        name, gtf_loc = line.split(",")
        genome = main(pathlib.Path(gtf_loc), 3000)
        genome[name]= genome
        
    


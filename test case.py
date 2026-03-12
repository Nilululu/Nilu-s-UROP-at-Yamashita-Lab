# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 15:41:16 2026

@author: nilof
"""

### tests the function "extract_genome_info()" 
###uses Guilt genome: the smallest genome among eaukaryotes 

from extract import extract_genome_info
from extract import get_genome_giant_introns

genomes= extract_genome_info("genomic.gtf")  #returns a dictionary of genomes


for genome in genomes.values():
     B= get_genome_giant_introns(genome, 3000)
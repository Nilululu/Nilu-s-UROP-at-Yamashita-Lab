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


genome_id, genes, giant_introns= main("genomic.gtf", 3000)

print(genome_id, giant_introns)

###gene_id "GUITHDRAFT_53115"; transcript_id "XM_005842043.1"

print(genes["GUITHDRAFT_53115"])

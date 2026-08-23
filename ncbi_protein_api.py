#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script is to download asnl files containing protein sequences and metadata using an API from NCBI

Created on Tue Aug 18 14:55:41 2026 by Nilu
"""

import pandas as pd
import requests
import logging 
import pathlib
import time

#mkaing a logger 
logging.basicConfig(filename = "ncbi_protein_logger.txt", level = logging.INFO, force = True)
logger = logging.getLogger(__name__)

#reading the table into a panda table
df = pd.read_csv("mammals_dystrophin.txt", sep = "\t", index_col = False)
df_sub = df[["species","protein_id"]]

#creating the folder that will store the protein asn1 files 
base_folder = pathlib.Path(pathlib.Path.cwd() / "ncbi_DMD_proteins")
base_folder.mkdir(exist_ok =  True)

#making sure we don't request more them 3 requests per second
start  =time.time()

#going through every species DMD protein_id in the table
for item in df_sub.itertuples(index = False):
    
    while time.time() - start < 0.4:
        time.sleep(0.1)
        
    name, protein_id = item
    
    #getting the file from NCBI
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={}&retmode=xml".format(protein_id)
    responce = requests.get(url) 
    
    # using species name and protein id for naming!
    file_name = name.replace(" ", "_") + "_" + protein_id + ".asn1"
    file = base_folder / file_name
    if responce.status_code == 200:
        with open(file, 'wb') as f:
            f.write(responce.content)
    else:
        logger.error("failed to downlaod the file with protein_id {}".format(protein_id))
        
    start = time.time()


    
    
            
    

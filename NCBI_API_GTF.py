# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 16:11:58 2026

@author: nilof
"""
import requests
import time
import pathlib
# some libraries later for extracting genome information out of each file


#maybe storing 10 accesion id from the table

with open ("ncbi_refseq-eukaryot.tsv", "r") as refseq_eukaryots:
    
    lines = refseq_eukaryots.readlines()

base_url = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GTF&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"
p = pathlib.Path.cwd()
base_dir = pathlib.Path(p/ "ncbi_data_directory")
base_dir.mkdir(exist_ok=True)  

start = time.time()
for i, line in enumerate(lines[1:]):
    if i >=10:
        break
        #find a way to skip the first line
    field= line.split("\t")
    id_= field[1]  # ge the id from the line
    name = field[3]
    while time.time() - start < 0.3:
        time.sleep(0.1)

     # ge the id from the line
    this_url = base_url.format(id_)
    response = requests.get(this_url)
    file_Path = base_dir /  '{}_genome.zip'.format(name)

    if response.status_code == 200:
        with open(file_Path, 'wb') as file:
            file.write(response.content)
        print('File downloaded successfully')
    else:
        print('Failed to download file')
    start = time.time()
        
# next steps: probably creating a function to unzip the files, go get the gtf file and add it to the genomes dictionary


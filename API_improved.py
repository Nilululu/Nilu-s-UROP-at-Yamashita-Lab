# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 15:39:01 2026

@author: nilof
"""

import requests
import time
import pathlib
import csv
from extract import main
# some libraries later for extracting genome information out of each file

step_2= 25  # is supposed to be 800
step = 5   #is supposed to be 10

genomes = dict()
# base url for downloading and base directory for storing
base_url = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GTF&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"
p = pathlib.Path.cwd()
base_dir = pathlib.Path(p/ "ncbi_data_directory")
base_dir.mkdir(exist_ok=True)  

# to make sure we don't request more then 5 urls per second we keep track of time
start = time.time()

directory_dict = dict()
with open ("ncbi_refseq-eukaryot.tsv", "r") as refseq_eukaryots:
    
    lines = refseq_eukaryots.readlines()

    for line_n, line in enumerate(lines[1:]):
        main_folder = line_n // step
        sub_folder = line_n % step
        top = line_n // step_2
        
        # makdir a folder with top_folder
        top_folder_path =  base_dir / 'top_{}'.format(str(top))
        top_folder_path.mkdir(exist_ok=True)
        
        # makdir a folder with main_fodler
        main_folder_path = top_folder_path / 'main_{}'.format(str(main_folder))
        main_folder_path.mkdir(exist_ok=True)
        
        # download the file in the main_folder
        field= line.split("\t")
        id_= field[1]  # ge the id from the line
        name = field[3]
        while time.time() - start < 0.3:
            time.sleep(0.1)
    
         # ge the id from the line
        this_url = base_url.format(id_)
        response = requests.get(this_url)
        file_path = main_folder_path /  '{}_genome.zip'.format(name)
    
        if response.status_code == 200:
            with open(file_path, 'wb') as file:
                file.write(response.content)
            print('File downloaded successfully')
        else:
            print('Failed to download file')
        start = time.time()
        
        dict_key= name + "_genome"
        dict_val= file_path
        directory_dict[dict_key]= dict_val
        print(top, main_folder, sub_folder)
        
        ### I am trying to extract info from each file as I go, but it still needs more work on
        
        # genome_file = list(file_path.rglob("genomic.gtf"))
        # extracted_info = main(genome_file[0], 3000)
        # genomes[name]= extracted_info
        
        
        
        if line_n > 20:  # for testing porpuses 
            break
    print(directory_dict)
    
    

    csv_file_name = "genomic_directory.csv"

    with open(csv_file_name, mode = 'w', newline= '') as file:
        writer= csv.writer(file)
        for key in directory_dict:
            writer.writerow([key, directory_dict[key]])  
            
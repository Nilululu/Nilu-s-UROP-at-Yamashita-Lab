# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 15:39:01 2026

@author: nilof
"""

import requests
import time
import pathlib
import csv
import zipfile

step_2 = 25  # is supposed to be 800
step = 5   #is supposed to be 10
directory_dict = dict() #to store genome names and genomic.gtf locatiosn
p = pathlib.Path.cwd() #defining the current folder
start = time.time()
id_set= set()


# base url for downloading 
base_url = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GTF&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"

#base directory for storing
base_dir = pathlib.Path(p/ "ncbi_data_directory")
base_dir.mkdir(exist_ok=True)   #creating the directory if not already created


with open ("ncbi_refseq-eukaryot.tsv", "r") as refseq_eukaryots:
    
    lines = refseq_eukaryots.readlines()

    for line_n, line in enumerate(lines[1:]):
        
        # download the file in the main_folder
        field= line.split("\t")
        id_= field[1]  # ge the id from the line
        if id_ in id_set:
            print("duplicated_id:", id_)
            continue

        id_set.add(id_)
        name = field[3]
        while time.time() - start < 0.3:
            time.sleep(0.1)
            
        main_folder = line_n // step
        sub_folder = line_n % step
        top = line_n // step_2
        
        # makdir a folder with top_folder
        top_folder_path =  base_dir / 'top_{}'.format(str(top))
        top_folder_path.mkdir(exist_ok=True)
        
        # makdir a folder with main_fodler
        main_folder_path = top_folder_path / 'main_{}'.format(str(main_folder))
        main_folder_path.mkdir(exist_ok=True)
        
    
        
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
        
        print(top, main_folder, sub_folder)
        
        zip_path = file_path
        extract_dir = file_path.parent / file_path.name.replace(".zip", "")
        
        
        with zipfile.ZipFile(zip_path, 'r') as z:
            z.extractall(extract_dir)
        
        
        # Search inside extracted folder for genomic.gtf
        genome_file = list(extract_dir.rglob("genomic.gtf"))
        print("Found GTF:", genome_file[0])
        
        # Handle potential missing GTF (should not accur though)
        if not genome_file:
            print(f"No genomic.gtf found for {name}, skipping")
        else:
            #save the location of the genomic file 
            dict_key = name + "_genome"
            dict_val = genome_file[0]
            directory_dict[dict_key]= dict_val
        if line_n > 100:  # for testing porpuses 
            break
    

    csv_file_name = "genomic_directory.csv"

    with open(csv_file_name, mode = 'w', newline= '') as file:
        writer= csv.writer(file)
        for key in directory_dict:
            writer.writerow([key, directory_dict[key]])  
            
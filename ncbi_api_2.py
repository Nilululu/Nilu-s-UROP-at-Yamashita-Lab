# -*- coding: utf-8 -*-
"""

The script uses an API to dowloand all ncbi eukaryotic annotated genomes ,excluding their FASTA files, from NCBI, 
The annotation can be both by refreq data base and individual submitters annotations 
The script keeps track of succssfully downloaded gtf files and stores their location in genomic_directory 

Created on July 16th by nilu
"""

import logging

logger = logging.getLogger(__name__)
logging.basicConfig(filename="error_download_log_2.txt", level=logging.WARNING, force = True)


import requests
import time
import pathlib
import csv
import zipfile
import shutil 




step_2 = 1000  # is supposed to be 800
step = 100   #is supposed to be 10

assert step_2 % step == 0

directory_dict = {} #to store genome names and genomic.gtf locatiosn
base_wkdir = pathlib.Path.cwd() 
start = time.time()
id_set = set()


# base url for downloading 
base_url = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{}/download?include_annotation_type=GENOME_GTF&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"

#base directory for storing
base_dir = pathlib.Path(base_wkdir/ "ncbi_data_directory_gca")
base_dir.mkdir(exist_ok = True)   #creating the directory if not already created


def download(base_url, id_, folder):
    
    """
    incorporates genome id into the Application Programming Interfase (API) url 
    and send a request to ncbi to download the zip folder in the given location

    Parameters
    base_url : str : NCBI API containing a request for different genome files that we need
    id_ : str : genome accession id 
    folder : path : path to the storage location 

    Returns
    -------
    file_path : path to the stored zip file

    
    """
    
    this_url = base_url.format(id_)
    response = requests.get(this_url)
    file_path = folder /  '{}.zip'.format(name)
    
    
        
    if response.status_code == 200:    # download the file in the main_folder
        with open(file_path, 'wb') as file:
            file.write(response.content)
        logger.info('File downloaded successfully')
    else:
        raise AssertionError
    return file_path
    


line_n = 0



with open ("ncbi_refseq-eukaryot.tsv", "r") as refseq_eukaryots:
    
    refseq_eukaryots.readline()
    for line in refseq_eukaryots:
        line = line.strip()
        if not line:
            continue

        field = line.split("\t")
        ids = field[1:3]
        id_ = None
        

        ids =  [ elem for elem in ids  if elem ]
        ids = sorted(ids)
        if len(ids) == 2:    #skipping genomes that have  two annotations (refseq + independant)
            continue
        
        id_ = ids[0]
        if id_.startswith("GCF"):   #skipping genomes that only have a refseq annotation
            continue 
        
        id_set.add(id_)

        name = field[3].replace(" ", "_") + "_" + id_
        while time.time() - start < 0.3:
            time.sleep(0.1)
            
        main_folder = line_n // step
        top = line_n // step_2

        # makdir a folder with top_folder
        top_folder_path =  base_dir / 'top_{}'.format(str(top))
        top_folder_path.mkdir(exist_ok=True)
        
        # makdir a folder with main_fodler
        main_folder_path = top_folder_path / 'main_{}'.format(str(main_folder))
        main_folder_path.mkdir(exist_ok=True)
        
        ### this is not working, can show current ones as suppressed and suppressed ones as current!
        # params = {"filters.assembly_version": "current"} # only latest, non-suppressed } 
        # params_response = requests.get( f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{id_}/dataset_report", params=params )
        
        # if params_response.status_code != 200:
            
        #     logger.warning(f"prevented downloading a suppressed GTF with name_id {name}")
        #     id_set.remove(id_)
        #     continue 
            
        
        try:
            file_path = download(base_url=base_url, id_=id_, folder=main_folder_path)
        except:
            logger.error(f"Failed to download zip file {id_}, {name}\n")
            id_set.remove(id_)
            continue
    
        # time after download used to  ward against API limit
        start = time.time()
        
        logger.info(top, main_folder)
        
        zip_path = file_path
        extract_dir = file_path.parent / file_path.name.replace(".zip", "")
        
        #unziping the downloaded file
        try:
            with zipfile.ZipFile(zip_path, 'r') as z:
                z.extractall(extract_dir)
        except zipfile.BadZipFile: 
            logger.error(f"Bad zipfile, {id_}, {name}\n")
            
            file_path.unlink()
            if extract_dir.is_dir():
                shutil.rmtree(extract_dir)
            id_set.remove(id_)
            continue 
        
        except:
            logger.error("enexpected error happened while unzipping the zip file")
            file_path.unlink()
            if extract_dir.is_dir():
                shutil.rmtree(extract_dir)
            id_set.remove(id_)
            continue

        file_path.unlink()
        
        # Search inside extracted folder for genomic.gtf
        genome_file = list(extract_dir.rglob("genomic.gtf"))

        
        # Handle potential missing GTF 
        if not genome_file:
            
            logger.error(f"No genomic.gtf found, {id_}, {name}\n")
            shutil.rmtree(extract_dir)
            id_set.remove(id_)
            continue
            

        logger.info("Found GTF:", genome_file[0])
        dict_key = genome_file[0]
        dict_val = name 
        directory_dict[dict_key]= dict_val
        line_n += 1
        
        # if line_n > 30:  # for testing purpuses 
        #    break        
       

csv_file_name = "genomic_directory_gca.csv"

with open(csv_file_name, mode = 'w', newline= '') as file:
    writer = csv.writer(file)
    for key in directory_dict:
        writer.writerow([key, directory_dict[key]])  

# save the id_set 
with open("id_set_gca.txt", "w") as f_id:
    
    for item in id_set:
        f_id.write(item)
        f_id.write("\n")
        
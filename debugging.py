# -*- coding: utf-8 -*-
"""
Created on Sun May 24 13:20:02 2026

@author: nilof
"""


def find_unique_GCF_from_path (filePath):
    """ a function that finds the unique words that starts with GCF_* in a text file and return them in a list/set"""
    
    gcfSet = set()
    with open (filePath, 'r') as file:
        
        file = file.readlines()
        for line in file:
            fields = line.split("/")
            for word in fields: 
                if word.startswith("GCF_"):
                    gcfSet.add(word)
    return gcfSet

def find_unique_GCF_column (filePath, n):
    """ a function that finds the unique words that starts with GCF_* in a text file and return them in a list/set"""
    
    gcfSet = set()
    with open (filePath, 'r') as file:
        
        for line in file:
            fields = line.split("\t")
            
            if (fields[n]).startswith("GCF_"):
                gcfSet.add(fields[n].strip())
            elif (fields[n]).startswith("GCA_"):
                continue
            #else:
                #print (fields[n])
  
    return gcfSet


def find_all_GCF (filePath):
    numgcf = 0
    gcf_list = []
    dashless_list = []
    with open (filePath, 'r') as file:
        
        for line in file:
            fields = line.split()
            
            for word in fields: 
                if word.startswith("GCF_"):
                    numgcf = numgcf + 1
                    gcf_list.append(word.strip())
                    dashless = (word.split("-"))[0]
                    dashless_list.append(dashless)
    print("num gcf apearance", numgcf)
    print("dashless set", len(set(dashless_list)))
    print("set", len(set(gcf_list)))
    
    return gcf_list, set(dashless_list)
    
    
    
    
    
def return_difference (setA, setB, setAname, setBname):
    
    A = len(setA)
    B = len(setB)
    
    print (setAname, A,"item")
    print(setBname, B, "item")
    
    print(f"items in {setAname} that are not in {setBname} include:{len(setA - setB)}items")
    print(setA - setB)
    
    print(f"items in {setBname} that are not in {setAname} include: {len(setB - setA)} item")
    
    print(setB - setA)
    print("\n"*2)

# # id_set = find_unique_GCF_column("id_set.txt", 0)    #should be replaced with id_set from cluster if obtained later
tsv = find_unique_GCF_column("ncbi_refseq-eukaryot.tsv", 1)
reports = find_unique_GCF_from_path("all_report_downloaded.txt")
genomic_directory = find_unique_GCF_from_path("genomic_directory.csv")
all_gtf_downloaded = find_unique_GCF_from_path("all_gtf_downloaded.txt")
id_set = find_unique_GCF_column("id_set.txt", 0)

print("tsv 2nd columns unique GCF", len (tsv))
print("genomic directory number", len(genomic_directory))
print("all_gtf_downloaded", len(all_gtf_downloaded))
print("id_set", len(id_set))
print("reports_downloaded", len(reports))

print("\n"*2)

return_difference(tsv, all_gtf_downloaded, "tsv", "downloaded_gtf")

return_difference(genomic_directory, all_gtf_downloaded, "genomic_directory", "all_gtf_downloaded")

return_difference(id_set, genomic_directory, "id_set", "genomic_directory")

return_difference(id_set, reports, "ide_set", "reports")

return_difference(reports, all_gtf_downloaded, "reports", "all_gtf")


    

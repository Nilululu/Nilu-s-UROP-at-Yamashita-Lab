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
        print("number of lines", len(file))
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
                gcfSet.add(fields[n])
            elif (fields[n]).startswith("GCA_"):
                continue
            else:
                print (fields[n])
  
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
                    gcf_list.append(word)
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
    
    print(f"items in {setAname} that are not in {setBname} include:")
    print(setA - setB)
    
    print(f"items in {setBname} that are not in {setAname} include:")
    
    print(setB - setA)

# # id_set = find_unique_GCF_column("id_set.txt", 0)    #should be replaced with id_set from cluster if obtained later
tsv_column_set = find_unique_GCF_column("ncbi_refseq-eukaryot.tsv", 1)
print("tsv 2nd columns unique GCF", len (tsv_column_set))

print("tsv comprehensive")
tsv_comprehensive_list, tsv_dashless_set = find_all_GCF("ncbi_refseq-eukaryot.tsv")

print("bash comprehensive")
tsv_bash_comprehensive_list, tsv_bash_dashless_set = find_all_GCF("bash_all.txt")

print("tsv bash unique")
tsv_bash_u_list, tsv_bash_u_dashless = find_all_GCF("bash_u.txt")

print("\nconfirmed that unique GCFs in bash take '-' into acount")


print("the difference between tsv second column unique variables and tsv comprehensive unique variables is about 115, why?")

return_difference(tsv_column_set, tsv_bash_dashless_set, "tsv 2nd column", "tsv comprehensive")

genomic_directory = find_unique_GCF_from_path("genomic_directory.csv")
all_gtf_downloaded = find_unique_GCF_from_path("all_gtf_downloaded_find.txt")

print("genomic directory number", len(genomic_directory))
print("all_gtf_downloaded", len(all_gtf_downloaded))


return_difference(tsv_column_set, all_gtf_downloaded, "tsv 2nd column", "downloaded_gtf")

return_difference(genomic_directory, all_gtf_downloaded, "genomic_directory", "all_gtf_downloaded")
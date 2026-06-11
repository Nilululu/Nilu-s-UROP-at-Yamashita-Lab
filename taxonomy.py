# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 14:35:27 2026

@author: nilof
"""
from pathlib import Path 

def generate_taxonomy_dict (nodes_file):
    
    """
    Parameters
    ----------
    nodes_file : dmp file containing taxIds. paren taxIds and rank
        
    Returns
    -------
    taxonomy dictionary : taxIds as keys, parent taxIds and ranks as values
    """
    
    taxonomy_dict = dict()
    with open (nodes_file, 'r') as nodes:
        lines = nodes.readlines()
        for line in lines:
            field = line.split("\t|")
            tax_id = int(field[0])
            parent_tax_id = int(field[1])
            rank = field[2].strip()
            taxonomy_dict[tax_id] = (parent_tax_id, rank)
    return taxonomy_dict  


def generate_tax_to_name(names_file):
    
    """
    Parameters
    ----------
    nodes_file : dmp file containing taxIds and names
        
    Returns
    -------
    tax_to-name : taxIds as keys, names as values
    """
    tax_to_name = dict()        
    with open(names_file, 'r') as names:
        lines = names.readlines() 
        for line in lines:
            field = line.split("\t|")
            tax_id = int(field[0])
            name = field[1].strip()
            tax_to_name[tax_id]= name
    return tax_to_name
       

       
def find_taxonomy(tax_id, taxonomy_dict, tax_to_name, famline_dict = {}):
    """
    

    Parameters
    ----------
    tax_id : int
        Taxonomy_id of a genome 
    taxonomy_dict : pyhton dictionary 
        a dictionary with all possible tax_ids as keys and parent_tax_id and rank as values
    tax_to_name : pyhton doctionary
        a dictionary with all possible tax_ids as keys and names as values
        famline_dict = dictionary that is initiated to store values of interest
    Returns
    -------
    famline_dictionary

    """
    
    try:    
        parent_id, rank = taxonomy_dict[tax_id]
        name = tax_to_name [tax_id]
        famline_dict[rank] = name
        
            
    except:
        print("Error: tax_id not found in the dictionary")
        return 
    
   

    if parent_id == 1:
        famline_dict["kingdom"] = "NA"
        return famline_dict
    if rank == "kingdom":
        return famline_dict
        
    else :

        return find_taxonomy(parent_id, taxonomy_dict, tax_to_name, famline_dict)
        


# nodes_file = Path(r"C:\Users\nilof\OneDrive\Documents\Python Scripts\new_taxdump\nodes.dmp")
# names_file = Path(r"C:\Users\nilof\OneDrive\Documents\Python Scripts\new_taxdump\names.dmp")

# taxonomy_dict = generate_taxonomy_dict(nodes_file)
# tax_to_name = generate_tax_to_name(names_file)


# print(find_taxonomy(51338, taxonomy_dict, tax_to_name))
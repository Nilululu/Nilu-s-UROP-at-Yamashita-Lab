# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 15:57:45 2026

@author: nilof
"""
import matplotlib.pyplot as plt
# a function to open FASTA files call for different information for each file

def give_feature_distribution(gene_file, feature):
    
    #reading the gene file and filtering for the desired feature
    gene = open(gene_file, "r") 
    feature_sizes=[]   #a list to store the feature size

    for line in gene:
        fields= line.strip().split("\t")
        if len(fields) < 3:
            continue 
        if fields[2]== feature:
            length= int(fields[4])-int(fields[3])  #end - start position
            feature_sizes.append(length)
    
    #print(intron_sizes)  #for debugging
            
    sorted_feature_sizes= sorted(feature_sizes)
    
    #plooting into a histogram
    plt.hist(sorted_feature_sizes, 1000)
    
    title = "Histogram for distribution of %s sizes" %feature
    plt.title(title)
    plt.xlabel("%s sizes"%feature)
    plt.ylabel("frequency")
    
    
give_feature_distribution("ACTA.txt", "intron")

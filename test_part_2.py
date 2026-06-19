#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 14:17:23 2026

connects the text file outputs of test_part_1 in a single text file and reads 
the data in that file to create various graphes
"""

from matplotlib import pyplot as plt
from plot_creator import get_style
from math import log10
import logging
import shutil

logging.basicConfig(filename="info_test2_logger.txt", level=logging.INFO)
logger = logging.getLogger(__name__)

# adding column names to the table
info_line = "#genome_id name kingdom tax_id total_sequence_length assembly_level assembly_type numChr num_scaffolds num_contigs scaffold_n50 contig_n50 gc_percent max_intron min_intron mean_intron median_intron sd_intron q_25 q_50 q_75 q_95 q_99 q_999 q_9999 q_99999"
info_line.replace(" ", "\t")

print(info_line)
text_file_list = ["table_0_6.txt", "table_6_12.txt", "table_18_24.txt", "table_12_18.txt", "table_24_30.txt"]
   
with open ("table.txt", "w") as destination:
    destination.write(info_line)
    destination.write("\n")
    for table in text_file_list:
        with open (table, 'r') as source:
            shutil.copyfileobj(source, destination)
    

#initializing some figues and dictionaries 

fig, ax = plt.subplots(1)
fig1, ax1 = plt.subplots(1)
ax1.set_title("genome size vs max intron (scale = kbp)")
ax.set_title("genome size vs max intron (scale = log bp)")
kingdom_count = dict()
type_count = dict ()
status_count = dict()
kingdom_mapping = dict()


#helper function
def add_counter(item, dictionary):
    """ if an item is not in dictionary the value of the added item will be 1, 
    else the count will increas by 1"""
    
    if item in dictionary:
        dictionary[item] +=1
    else:
        dictionary[item] = 1 


#reading from the table text file
with open ("table.txt", 'r') as table_file:
    
    for line in table_file:
        if line.startswith("#"):  #skipping the info line
            continue 
        
        try:
            fields = line.strip().split("\t")
            assembly_status = fields[5]
            kingdom = fields[2]
            assembly_type = fields[6]
            size = int(fields[4])
            max_intron = int(fields[13])
            
            add_counter(assembly_status, status_count)
            add_counter(kingdom, kingdom_count)
            add_counter(assembly_type, type_count)
            
            # getting a random color for each kingdom in the graph
            get_style(kingdom, kingdom_mapping)
            style = kingdom_mapping[kingdom]
            ax1.scatter((size*0.001), (max_intron*0.001), color = style['color'])
            
            #ax1.legend()  # need to learn how to make legends, for now logging the legends
            
            ax.scatter(log10(size), log10(max_intron), color = style['color'])
           
            
        except:
            logger.error(f"failed to etract information for {line}")


# ax.legend(labels = kingdom_mapping, bbox_to_anchor= (1.05, 1))

handle = []
handle_1 = []
for item in kingdom_mapping:
    handle.append(ax.scatter([], [], color = kingdom_mapping[item]["color"], label = item))
    handle_1.append(ax1.scatter([], [], color = kingdom_mapping[item]["color"], label = item))

ax.legend(handles = handle)
ax1.legend(handles = handle_1)

#kingdom distribution of the data
fig2, ax2 = plt.subplots(1)
ax2.pie(kingdom_count.values(), labels = kingdom_count.keys())

#assembly type distribution of the data
fig3, ax3 = plt.subplots(1)
ax3.pie(type_count.values())
ax3.legend( loc='upper left', labels = type_count.keys(), bbox_to_anchor=(1.05, 1))
fig3.tight_layout()

#assembly status distribution 
fig4, ax4 = plt.subplots(1)
ax4.pie(status_count.values())

# Place the legend outside the chart (top-right corner area)
ax4.legend(loc='upper left', labels= status_count.keys(), bbox_to_anchor=(1.05, 1))

# Prevent the outside legend from getting cut off in the saved/rendered image
fig4.tight_layout()
print(kingdom_mapping)


#logging the dictionaries created
logger.info(kingdom_count)
logger.info(kingdom_mapping)
logger.info(status_count)
logger.info(type_count)

#saving the generated figures
fig.savefig("genome_size_vs_max_intron_log")
fig1.savefig("genome_size_vs_max_intron_kbp")
fig2.savefig("kingdom_distribution")
fig3.savefig("type_distribution") 
fig4.savefig("status_distributio")



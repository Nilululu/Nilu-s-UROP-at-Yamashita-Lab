#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 14:17:23 2026

connects the text file outputs of test_part_1 in a single text file and reads 
the data in that file to create various graphes
"""

from matplotlib import pyplot as plt
import random
from math import log10
import logging
import numpy as np


##### logging file created
logging.basicConfig(filename="info_test2_logger.txt", level=logging.INFO)
logger = logging.getLogger(__name__)


##### HELPER FUNCTIONS

def get_style (kingdom, _mapping = dict(), _used = set()):
    
    """
    Summery: tasigns a random color to every kingdom and stores it in _mapping
    
    Parameters:
        kingdom: str
        _mapping: python dictionary
        _used: set
    
    retruns
    _mapping (updated)
    
    """
    colors = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e',
     '#e6ab02','#a6761d','#666666']
    
    #returning existing info
    if kingdom in _mapping:
        return _mapping
    
    #looping until we get a color that has not been used before
    while True:
        color = random.choice(colors)
        if color in _used:
            continue
        
        _used.add(color)
        
        break 
        
    _mapping[kingdom] = {"color" : color}
    
    return _mapping
        
    
def add_counter(item, dictionary):
    """ if an item is not in dictionary the value of the added item will be 1, 
    else the count will increas by 1"""
    
    if item in dictionary:
        dictionary[item] +=1
    else:
        dictionary[item] = 1 
        



##### Initializing some figues for scatter plots and dictionaries 
fig, ax = plt.subplots(1)
ax.set_title("genome size vs max intron (scale = log bp)")
ax.set_xlabel("genome size (log10)")
ax.set_ylabel("max intron (log10)")

fig1, ax1 = plt.subplots(1)
ax1.set_title("genome size vs max intron (scale = kbp)")
ax1.set_xlabel("genome size (kbp)")
ax1.set_ylabel("max intron (kbp)")

kingdom_count = dict()
type_count = dict ()
status_count = dict()
kingdom_mapping = dict()

##### Reading data from a text file

with open ("result_table.txt", 'r') as table_file:
    
    boxplots_data = np.zeros((len(table_file.readlines()), 13))
    
    i = 0 #used to write data in numpy arrays
    for line in table_file:
        
        if line.startswith("#"):  #skipping the info line
            continue
        
        try:
            fields = line.strip().split("\t")
            
            #getting the information needed for graphes
            assembly_status = fields[5]
            kingdom = fields[2]
            assembly_type = fields[6]
            size = int(fields[4])
            max_intron = int(fields[13])
            
            
            # will contain max_intron, min_intron, mean_intron,median_intron, 
            #sd_intron, q_25, q_50, q_75, q_95, q_99, q_999, q_9999, q_99999
            boxplots_data[i, 0:13] = fields[13:26]
            
            #updating counters 
            add_counter(assembly_status, status_count)
            add_counter(kingdom, kingdom_count)
            add_counter(assembly_type, type_count)
            
            # getting a random color for each kingdom in the graph
            get_style(kingdom, kingdom_mapping)
            style = kingdom_mapping[kingdom]
            
            #plotting genome size vs max intron scale kbp
            ax1.scatter((size*0.001), (max_intron*0.001), color = style['color'])
            
            #plotting the same thing but with log 10 scale
            ax.scatter(log10(size), log10(max_intron), color = style['color'])
           
            i +=1
        except:
            logger.error(f"failed to extract information for {line}")
        
        # print(boxplots_data)
    #making color legends for fig, and fig1
    handle = []
    handle_1 = []
    for item in kingdom_mapping:
        handle.append(ax.scatter([], [], color = kingdom_mapping[item]["color"], label = item))
        handle_1.append(ax1.scatter([], [], color = kingdom_mapping[item]["color"], label = item))
    
    ax.legend(handles = handle)
    ax1.legend(handles = handle_1)
    


    ##### MAKING MORE PLOTS


    #pie chart for kingdom distribution in the data
    fig2, ax2 = plt.subplots(1)
    ax2.pie(kingdom_count.values(), labels = kingdom_count.keys())
    ax2.set_title("kingdom distribution")
    
    #assembly type distribution of the data
    fig3, ax3 = plt.subplots(1)
    ax3.pie(type_count.values())
    ax3.legend( loc='upper left', labels = type_count.keys(), bbox_to_anchor=(1.05, 1))
    ax3.set_title("assembly type distribution")
    fig3.tight_layout()
    
    #assembly status distribution 
    fig4, ax4 = plt.subplots(1)
    ax4.pie(status_count.values())
    ax4.legend(loc='upper left', labels= status_count.keys(), bbox_to_anchor=(1.05, 1))
    ax4.set_title("assembly status")
    fig4.tight_layout()
    
    #box plots for different quantiles + max and min with outliers
    fig5, ax5 = plt.subplots()
    ax5.boxplot((boxplots_data[:, [1,5,6,7,8,9,10,11,12,0]])*0.001)
    ax5.set_xticklabels(["min", "25", "50", "75", "95", "99", "999", "9999", "99999", "max"])
    ax5.set_xlabel("intron quantiles")
    ax5.set_title("length distribution of introns in each intron quantile across genomes")
    ax5.set_ylabel("Length (kbp)")
    
    #trying 10g10 scale
    fig6, ax6 = plt.subplots()
    ax6.boxplot(np.log10(boxplots_data[:, [1,5,6,7,8,9,10,11,12,0]]))
    ax6.set_xticklabels(["min", "25", "50", "75", "95", "99", "999", "9999", "99999", "max"])
    ax6.set_xlabel("intron quantiles")
    ax6.set_title("length distribution of introns in each intron quantile across genomes")
    ax6.set_ylabel("Length (log10)")
    
    #box plot for different quantiles + max and min without outliers
    fig7, ax7 = plt.subplots()
    ax7.boxplot((boxplots_data[:, [1,5,6,7,8,9,10,11,12,0]])*0.001, showfliers = False)  #from q_25 to q_99999
    ax7.set_xticklabels(["min", "25", "50", "75", "95", "99", "999", "9999", "99999", "max"])
    ax7.set_xlabel("intron quantiles")
    ax7.set_title("length distribution of introns in each intron quantile across genomes without outliers")
    ax7.set_ylabel("Length (kbp)")
   
    #trying log10 scale
    fig8, ax8 = plt.subplots()
    ax8.boxplot(np.log10(boxplots_data[:, [1,5,6,7,8,9,10,11,12,0]]), showfliers = False)  #from q_25 to q_99999
    ax8.set_xticklabels(["min", "25", "50", "75", "95", "99", "999", "9999", "99999", "max"])
    ax8.set_xlabel("intron quantiles")
    ax8.set_title("length distribution of introns in each intron quantile across genomes without outliers")
    ax8.set_ylabel("Length (log10)")
    
    ##### Logging the counter dictionaries created
    logger.info(kingdom_count)
    logger.info(kingdom_mapping)
    logger.info(status_count)
    logger.info(type_count)
    
    ##### Saving the generated figures
    fig.savefig("genome_size_vs_max_intron_log.png")
    fig1.savefig("genome_size_vs_max_intron_kbp.png")
    fig2.savefig("kingdom_distribution.png")
    fig3.savefig("type_distribution.png") 
    fig4.savefig("status_distribution.png")
    fig5.savefig("box_plots_ouliers_kbp.png")
    fig6.savefig("box_plots_outliers_log10.png")
    fig7.savefig("box_plots_no_oulier_kbp.png")
    fig8.savefig("box_plots_no_outlier_log10.png")




        






















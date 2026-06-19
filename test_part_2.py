#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 14:17:23 2026

@author: nilofarghafory
"""

from matplotlib import pyplot as plt
from plot_creator import get_style
from math import log10
import logging
import shutil

logging.basicConfig(filename="info_test2_logger.txt", level=logging.INFO)
logger = logging.getLogger(__name__)


# a function to join all the taxt files

#list of text files to join in order


text_file_list = ["table_0_6.txt", "table_6_12.txt", "table_18_24.txt", "table_12_18.txt", "table_24_30.txt"]

with open ("table.txt", "wb") as destination:
    
    for table in text_file_list:
        with open (table, 'rb') as source:
            shutil.copyfileobj(source, destination)
    


fig, ax = plt.subplots(1)
fig1, ax1 = plt.subplots(1)
ax1.set_title("genome size vs max intron (scale = kbp)")
ax.set_title("genome size vs max intron (scale = log bp)")
kingdom_count = dict()
type_count = dict ()
status_count = dict()
kingdom_mapping = dict()

with open ("table.txt", 'r') as table_file:
    
    for line in table_file:
        
        try:
            fields = line.strip().split("\t")
            status = fields[5]
            kingdom = fields[2]
            g_type = fields[6]
            size = int(fields[4])
            max_intron = int(fields[13])
            
            if status not in status_count:
                status_count[status] =1
            else:
                status_count[status] += 1
            
            if kingdom not in kingdom_count:
                kingdom_count[kingdom] = 1
            else:
                kingdom_count[kingdom] += 1
            
            if g_type not in type_count:
                type_count[g_type] = 1
            else:
                type_count[g_type] += 1
             
            get_style(kingdom, kingdom_mapping)
            style = kingdom_mapping[kingdom]
            ax1.scatter((size*0.001), (max_intron*0.001), color = style['color'])
            
            #ax1.legend()  # need to learn how to make legends, for now logging the legends
            
            ax.scatter(log10(size), log10(max_intron), color = style['color'])
            #ax.legend()
            
        except:
            logger.error(f"failed to etract information for {line}")

ax.legend(labels = kingdom_mapping, bbox_to_anchor= (1.05, 1))
fig2, ax2 = plt.subplots(1)
ax2.pie(kingdom_count.values(), labels = kingdom_count.keys())

fig3, ax3 = plt.subplots(1)
ax3.pie(type_count.values())

# Place the legend outside the chart (top-right corner area)
ax3.legend( loc='upper left', labels = type_count.keys(), bbox_to_anchor=(1.05, 1))
fig3.tight_layout()

fig4, ax4 = plt.subplots(1)
ax4.pie(status_count.values())

# Place the legend outside the chart (top-right corner area)
ax4.legend(loc='upper left', labels= status_count.keys(), bbox_to_anchor=(1.05, 1))

# Prevent the outside legend from getting cut off in the saved/rendered image
fig4.tight_layout()
print(kingdom_mapping)

logger.info(kingdom_count)
logger.info(kingdom_mapping)
logger.info(status_count)
logger.info(type_count)

fig.savefig("genome_size_vs_max_intron_log")
fig1.savefig("genome_size_vs_max_intron_kbp")
fig2.savefig("kingdom_distribution")
fig3.savefig("type_distribution") 
fig4.savefig("status_distributio")



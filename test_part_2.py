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

logging.basicConfig(filename="info_test2_logger.txt", level=logging.INFO)
logger = logging.getLogger(__name__)


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
            status = fields[2]
            kingdom = fields[8]
            g_type = fields[3]
            size = int(fields[4])
            max_intron = int(fields[7])
            
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

fig2, ax2 = plt.subplots(1)
ax2.pie(kingdom_count.values(), labels = kingdom_count.keys())

fig3, ax3 = plt.subplots(1)
ax3.pie(type_count.values(), labels = type_count.keys())

fig4, ax4 = plt.subplots(1)
ax4.pie(status_count.values(), labels = status_count.keys())

logger.info(kingdom_count)
logger.info(kingdom_mapping)
logger.info(status_count)
logger.info(type_count)

fig.savefig("genome_size_vs_max_intron_log")
fig1.savefig("genome_size_vs_max_intron_kbp")
fig2.savefig("kingdom_distribution")
fig3.savefig("type_distribution") 
fig4.savefig("status_distributio")
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 22:31:56 2017

This script takes a gene set and a ranked gene list and calculates
the ES for the gene set.

@author: thasegawa
"""

import pandas as pd

# Specify filenames
ranks_fname = 'P53_hgu95av2_ES.txt'
probe2gene_fname = 'P53_hgu95av2_genes_v2.txt'
gset_fname = 'chr3q29.gmt'

# Read probe rank file
ranks = pd.read_table(ranks_fname)

# Read probe 2 gene converter
probe2gene = {}
with open(probe2gene_fname) as f:
    lines = f.readlines()
    for line in lines:
        splits = line.strip('\n').split('\t')
        probe = splits[0]
        
        if len(splits) > 1:
            genes = splits[1:]
            
        probe2gene[probe] = genes

# Read geneset
gset = []
with open(gset_fname) as f:
    gset = f.read().split('\t')

# Cut off header from gene set file
gset = gset[2:]

# Go through each probe and calculate running sum
rsum = 0
count = 0
rsum_list = []
for [probe, score] in zip(ranks['NAME'], ranks['ES']):
    genes = probe2gene[probe]
    
    # Look for the gene in the gene set
    trigger = False
    for gene in genes:
        if gene in gset:
            trigger = True
            break
    
    # Increase the running-sum if the gene is in the gene set. If not,
    # decrease the running sum        
    if trigger:
        count += 1
        rsum += score
        print(gene)
    else:
        rsum -= score
        
    rsum_list.append(rsum)
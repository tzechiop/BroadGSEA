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
N = ranks.shape[0]

# Specify parameters
p = 1

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
N_H = len(gset)

# Find all genes in geneset to calculate N_R
N_R = 0
hitcount = 0
for [probe, rval] in zip(ranks['NAME'], ranks['R']):
    genes = probe2gene[probe]
    for gene in genes:
        if gene in gset:
            N_R += abs(rval)**p
            hitcount += 1
            break
        
print('Hitcount = {0}'.format(hitcount))
print('N_R = {0:.3f}'.format(N_R))

# Go through each probe and calculate running sum
count = 0
rsum = 0
pval = 0
rsum_list = []
pval_list = []
rank_list = []
hitgene_list = []
probe_list = []
hitrval_list = []
for index, (probe, rval) in enumerate(zip(ranks['NAME'], ranks['R'])):
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
        cur_pval = abs(rval)**p/N_R
        rank_list.append(index + 1)
        hitgene_list.append(gene)
        probe_list.append(probe)
        hitrval_list.append(rval)
    else:
        cur_pval = -1/(N - N_H)

    # Calculate running sum        
    rsum += cur_pval

    pval_list.append(cur_pval)        
    rsum_list.append(rsum)
    
ranks['Pval'] = pval_list
ranks['Rsum'] = rsum_list
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 22:53:48 2017

@author: thasegawa
"""

import csv
import mygene
import urllib

# Read expression file
fname_expr = 'P53_hgu95av2.gct'
expr = []
with open(fname_expr) as f:
    for line in csv.reader(f, dialect = "excel-tab"):
        expr.append(line)
        
# Cut off header from expression file
expr = expr[3:]

# Read gene set file
fname_gset = 'chr7p15.gmt'
gset = []
with open(fname_gset) as f:
    gset = f.read().split('\t')

# Cut off header from gene set file
gset = gset[2:]

# Iterate through probe sets and find corresponding gene
probes = {}
extraprobes = {}
mg = mygene.MyGeneInfo()
for index, line in enumerate(expr[12579:]):
    probe = line[0]
    query = probe.split('/')[0]
    result = mg.query(query)
    try:
        hits = result['hits']
        genes = []
        for hit in hits:
            genes.append(hit['symbol'])
        probes[probe] = genes
    except Exception as e:
        print('Error for probe #{0}: {1}'.format(index + 1, e))
    
    if ((index + 1) % 1000 == 0) and (index > 0):
        print('{0} probes analyzed!'.format(index + 1))
    
            
            
    
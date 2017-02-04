# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 22:53:48 2017

This script converts the probes contained within expression data to
genes using the mygene library. The script outputs a file that lists
each probes and its associated gene(s).

http://mygene.info/

@author: thasegawa
"""

import csv
import mygene

# Parameters
fname_expr = 'P53_hgu95av2.gct'
outfname = 'P53_hgu95av2_genes_v2.txt'

# Read expression file
expr = []
with open(fname_expr) as f:
    for line in csv.reader(f, dialect = "excel-tab"):
        expr.append(line)

# Iterate through probe sets and find corresponding gene
probes = {}
mg = mygene.MyGeneInfo()
for index, line in enumerate(expr    ):
    probe = line[0]
    query = probe.split('/')[-1]
    result = mg.query(query)
    try:
        hits = result['hits']
        genes = []
        for hit in hits:
            genes.append(hit['symbol'])
        probes[probe] = genes
    except Exception as e:
        print('Error for probe #{0}:'.format(index + 1))
        print(e)
    
    if ((index + 1) % 1000 == 0) and (index > 0):
        print('{0} probes analyzed!'.format(index + 1))
    
# Output probe-gene file
with open(outfname, 'w') as f:
    for probe, genes in probes.items():
        f.write('{0}\t{1}\n'.format(probe, '\t'.join(genes)))

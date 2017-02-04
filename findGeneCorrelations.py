# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 22:05:51 2017

This script goes through the genes in a probe and calculates their expression score.

@author: thasegawa
"""

import pandas as pd
import numpy as np

# Function to calcuate the ranking score see page 43 in GSEA User Guide
def calcScore(row):
    vals = {pheno: [] for pheno in uniqpheno}
    for expr, pheno in zip(row[2:], phenotypes):
        vals[pheno].append(expr)
    u_0 = np.mean(vals[uniqpheno[0]])
    u_1 = np.mean(vals[uniqpheno[1]])
    s_0 = np.std(vals[uniqpheno[0]])
    s_1 = np.std(vals[uniqpheno[1]])
        
    ES = (u_0 - u_1)/(s_0 + s_1)
    
    return ES
    
# Specify filenames and phenotype labels
fname_pheno = 'P53.cls'
fname_expr = 'P53_hgu95av2.gct'
outfname = 'P53_hgu95av2_ES.txt'
uniqpheno = ['WT', 'MUT']
outcols = ['NAME', 'Score']

# Read in phenotypes of each sample
with open(fname_pheno) as f:
    lines = f.readlines()
phenotypes = lines[2].strip().split(' ')

# Read in expression data
expr = pd.read_table(fname_expr, header = 2)

# Calculate expression score
expr['Score'] = expr.apply(calcScore, axis = 1)

# Sort by expression score and output
expr = expr.sort_values('ES', ascending = False)
expr[outcols].to_csv(outfname, sep = '\t', index = False)
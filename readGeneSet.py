# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 23:56:47 2017

@author: thasegawa
"""

# Cut off header from expression file
expr = expr[3:]

# Read gene set file
fname_gset = 'chr7p15.gmt'
gset = []
with open(fname_gset) as f:
    gset = f.read().split('\t')

# Cut off header from gene set file
gset = gset[2:]
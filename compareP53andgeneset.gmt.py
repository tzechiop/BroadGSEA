# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 22:53:48 2017

@author: thasegawa
"""

import urllib

def uniprot_mapping(fromtype, totype, identifier):
    base = 'http://www.uniprot.org'
    tool = 'mapping'
    params = {'from':fromtype,
                'to':totype,
                'format':'tab',
                'query':identifier,
    }
    data = urllib.urlencode(params)
    url = base+'/'+tool+'?'+data
    response = urllib.request(url)
    return response.read()
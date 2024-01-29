#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 01:30:19 2024

@author: samiha
"""
import pandas as pd
import numpy as np

def rounding(n):
    if (type(n)==str) or (n==None) or (np.isnan(n)):
        return str(n)
    elif n==0:
        return str(0)
    elif ((abs(n)>=1e-2) and (abs(n)<=1e4)):
        sgn = '-' if n<0 else ''
        num = format(abs(n) - int(abs(n)),'f')
        if int(num[2:])<1:
            d = str(abs(n))
            return sgn + d[:-2]
        else:
            return str(np.round(n,2))
    else:
        return '{:.0e}'.format(n)

def read_rounding(n):
    if (n=='nan') or (n=='None'):
        return np.nan
    elif (type(n)==str):
        try:
            convert2float = float(n)
            return convert2float
        except:
            return str(n)

# write or display
#aa = rdf.applymap(lambda x: rounding(x))
# read for calc.
#bb = aa.applymap(lambda x: read_rounding(x))        

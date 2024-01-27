#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 19:54:16 2024

@author: samiha
"""

import pandas as pd
import numpy as np

def geos_chem_species(nmog):    
    # S07 model species of unassigned compounds
    ll_s07=['ALK4','ALK5','BALD','CRES', 'OLE1', 'OLE2', 'RCHO']
    # S18 species
    ll_s18=['FURNS', 'XYNL', 'SVPHE','NAPS']
    # Corresponding GEOOS_chem model species of S07
    g_s07=['ALK4','ALK4','BALD','CSL', 'PRPE', 'PRPE', 'RCHO']
    # Corresponding GEOOS_chem model species of S18B
    g_s18=['FURA', 'MCT','MCT','NAPS']
    
    # explicit spc
    spc=['C2H2','CH2O']
    com=['Ethyne', 'Formaldehyde']
    
    adf=pd.DataFrame()
    adf['MOZRT_sp']=ll2
    adf['S18B_sp']=ll
    
    # Assign MOZT1 model species
    for aa in range(len(adf)):
        s=adf['S18B_sp'].iloc[aa]
        m=adf['MOZRT_sp'].iloc[aa]
        ind=nmog[nmog['S18B']==s][nmog['MOZT1'].isnull()].index
        if len(ind)!=0:
            for ii in ind:
                nmog.loc[ii,'MOZT1']=m
    return nmog


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
    ll_s07=['HCHO','ALK4','ALK5','BALD','CRES', 'OLE1', 'OLE2', \
            'RCHO', 'ARO1', 'ARO2','TERP','MVK','MACR','ACYE','PACD',\
               'IPRD', 'BACL', 'PRD2', 'ALK3', 'ALK1', 'ALK2', 'RNO3', 'ROOH', 'MEK']
    # S18 species
    ll_s18=['FURNS','PHEN','NAPS']
    # Corresponding GEOOS_chem model species of S07
    g_s07=['CH2O','ALK4','ALK4','BALD','CSL', 'PRPE', \
           'PRPE', 'RCHO', 'TOLU', 'XYLE', 'MTPO','MVK','MACR','ALK4','ALK4',\
               'MVK', 'MGLY','MEK', 'ALK4', 'C2H6', 'C3H8', 'R4N2', 'MP', 'MEK']
    # Corresponding GEOOS_chem model species of S18B
    g_s18=['FURA','PHEN','NAPS']
    
    adf=pd.DataFrame()
    adf['geos_chem_sp']=g_s07
    adf['S07_sp']=ll_s07
    
    adf2=pd.DataFrame()
    adf2['geos_chem_sp']=g_s18
    adf2['S18B_sp']=ll_s18
        
    # Assign GEOS_chem model species
    for aa in range(len(adf)):
        s=adf['S07_sp'].iloc[aa]
        g=adf['geos_chem_sp'].iloc[aa]
        ind=nmog[nmog['S07']==s][nmog['GEOS_chem'].isnull()].index
        if len(ind)!=0:
            for ii in ind:
                nmog.loc[ii,'GEOS_chem']=g
                
    # Assign GEOS_chem model species
    for aa in range(len(adf2)):
        s=adf2['S18B_sp'].iloc[aa]
        g=adf2['geos_chem_sp'].iloc[aa]
        ind=nmog[nmog['S18B']==s][nmog['GEOS_chem'].isnull()].index
        if len(ind)!=0:
            for ii in ind:
                nmog.loc[ii,'GEOS_chem']=g
                
    return nmog

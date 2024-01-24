#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 15:32:46 2024

@author: samiha
"""

import pandas as pd
import numpy as np

def mozart_species(nmog):
    # Unassigned MOZT1 and mapping S18B species to MOZT1 species
    nmog['MOZT1'].unique()
    nmog['S18B'][nmog['MOZT1'].isnull()].unique()
    
    # S18 model species of unassigned compounds
    ll=['HCHO','KET2','OTH3','AMINS','ACYLS','OTH4','NROG',
        'OLEP','ROOH','OLE1', 'OLE2','FURNS','RCOOH','RCHO',
        'R1NO3','LVKS','OLEC','ARO1','OLED','BUT13','OLEA2',
        'OLEA1','AFG1','STYRS','PHEN','OLE4','ALK5','MEK','CRES',
        'ARO2','OTH2','NAPS','DLIMO','SVPHE','MVK','XYNL', 'SESQ']
    
    # Corresponding MOZT1 model species
    ll2=['CH2O','MEK','BIGALK','C3H6','BIGENE','BIGENE','NROG',
         'BIGENE','XYLENES', 'BIGENE', 'BIGENE','XYLENES','TOLUENE',
         'BIGALK','ALKNIT','MEK','BIGENE','XYLENES','BIGENE','BIGENE','MEK',
         'TOLUENE','MEK','BIGENE','PHENOL','BIGENE','BIGALK','MEK','CRESOL',
         'XYLENES','NROG','XYLENES','LIMON','PHENOL','MVK','CRESOL','BCARY']
    
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


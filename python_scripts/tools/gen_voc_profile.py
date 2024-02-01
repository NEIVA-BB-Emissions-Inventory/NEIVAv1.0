#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 11:53:38 2023

@author: samiha
"""

import pubchempy as pcp
import pandas as pd
import numpy as np
from NEIVA.python_scripts.connect_with_mysql import*
from NEIVA.python_scripts.tools.assign_mozart_species import mozart_species
from NEIVA.python_scripts.tools.assign_geos_chem_species import geos_chem_species
from NEIVA.python_scripts.tools.join_ef_property_table import *
from NEIVA.python_scripts.tools.number_format_function import *
from sqlalchemy import text
import matplotlib.pyplot as plt

            
def voc_profile(dd, chem, ft):
    output_db=connect_db('neiva_output_db')
    bk_db=connect_db('backend_db')
    
    nmog=join_ef_property(dd)
    
    # Set EF column based on the input parameter 'fire type'
    efcol='AVG_'+ft.replace(' ','_')
    nmog['ef']=nmog[efcol]
    nmog=nmog[nmog['ef'].notna()].reset_index(drop=True)

    nmog_final = lump_com_with_speciation (nmog, chem)
    nmog_final = distribute_unk_ef (dd, efcol, nmog_final)
    
    nmog_final['mole']=nmog_final['ef']/nmog_final['mm']
    # The list of unique model species
    uu=list(nmog_final[chem][nmog_final[chem].notna()].unique())
    
    # Creating the model species profile dataframe
    prdf=pd.DataFrame()
    for i in range(len(uu)):
        prdf.loc[i,chem]=uu[i]        
        totef=nmog_final['ef'][nmog_final[chem]==uu[i]].sum()
        mean_mm=nmog_final['mm'][nmog_final[chem]==uu[i]].mean()
        prdf.loc[i,'ef']=totef
        prdf.loc[i,'mm']=mean_mm
                
    prdf['mole']=prdf['ef']/prdf['mm']
    
    totmole=prdf['mole'].sum()
    prdf['mole_fraction']=prdf['mole']/totmole
    prdf=prdf.sort_values(by='mole_fraction', ascending=False)
    prdf=prdf.reset_index(drop=True)
    # prdf['mole']=round(prdf['mole'],4)    
    # prdf['mole_fraction']=round(prdf['mole_fraction'],4)
    prdf = prdf.applymap(lambda x: rounding(x))
    return prdf

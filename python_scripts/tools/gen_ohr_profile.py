#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 21:00:30 2023

@author: samiha
"""

import pubchempy as pcp
import pandas as pd
import numpy as np
from NEIVA.python_scripts.connect_with_mysql import*
from NEIVA.python_scripts.tools.assign_mozart_species import mozart_species
from NEIVA.python_scripts.tools.join_ef_property_table import *
from NEIVA.python_scripts.tools.assign_geos_chem_species import geos_chem_species
from NEIVA.python_scripts.tools.number_format_function import *
from sqlalchemy import text
import matplotlib.pyplot as plt


def calc_OHR (dd,chem,ft, tot_voc):
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
    
    nmog_final['mole_frac']=nmog_final['mole']/nmog_final['mole'].sum()
    nmog_final['conc']=nmog_final['mole_frac']*tot_voc # tot_voc in ppb
    nmog_final['ohr']=nmog_final['conc']*2.5e10*nmog_final['kOH']      

    uu=list(nmog_final[chem][nmog_final[chem].notna()].unique())
    # Creating the model species profile dataframe
    ohrdf=pd.DataFrame()
    for i in range(len(uu)):
        ohrdf.loc[i,chem]=uu[i]
        ohrdf.loc[i,'OHR']=nmog_final['ohr'][nmog_final[chem]==uu[i]].sum()
    
    tot_ohr=ohrdf['OHR'].sum()
    ohrdf['OHR_frac']=ohrdf['OHR']/tot_ohr
    
    ohrdf=ohrdf.sort_values(by='OHR', ascending=False)
    ohrdf=ohrdf.reset_index(drop=True)
    
    ohrdf=ohrdf.applymap(lambda x: rounding(x))

    return ohrdf

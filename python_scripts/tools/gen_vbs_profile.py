#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 19:30:13 2023

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


def calc_VBS (dd, ft):
    output_db=connect_db('neiva_output_db')
    bk_db=connect_db('backend_db')
    
    nmog=join_ef_property(dd)
    
    # Set EF column based on the input parameter 'fire type'
    efcol='AVG_'+ft.replace(' ','_')
    nmog['ef']=nmog[efcol]
    nmog=nmog[nmog['ef'].notna()].reset_index(drop=True)
    
    chem='S07'
    nmog_final = lump_com_with_speciation (nmog, chem)
    nmog_final = distribute_unk_ef (dd, efcol, nmog_final)
            
    nmog_final['bin']=pd.cut(nmog_final['cstar'], bins=range(0,14,1))
    aa=nmog_final.groupby('bin')['ef'].sum()/nmog_final['ef'].sum()
    
    fdf=pd.DataFrame()
    fdf['bin']=aa.index
    fdf['ef/sum_ef']=aa
    fdf['ef/sum_ef']=fdf['ef/sum_ef'].map(lambda x:rounding(x))
    return fdf





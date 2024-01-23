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

    nmog_final = lump_com_with_speciation (nmog, chem)
    nmog_final = distribute_unk_ef (dd, efcol, nmog_final)
    
    nmog_final['mole']=nmog_final['ef']/nmog_final['mm']
    
    
    ft=ft.replace(' ','_')
    efcol='AVG_'+ft
    nmog['ef']=nmog[efcol]
    
    nmog=nmog[nmog['ef'].notna()].reset_index(drop=True)
    
    nmog_final = lump_com_with_speciation (nmog, chem)
    nmog_final = distribute_unk_ef (dd, efcol, nmog_final)
    
    nmog_final['bin']=pd.cut(nmog_final['cstar'], bins=range(0,14,1))
    aa=nmog_final.groupby('bin')['ef'].sum()/nmog_final['ef'].sum()
    
    fdf=pd.DataFrame()
    fdf['bin']=aa.index
    fdf['ef/sum_ef']=aa

    return fdf





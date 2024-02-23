#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 11:45:35 2024

@author: samiha
"""

import pandas as pd
import numpy as np
import pubchempy as pcp
import matplotlib.pyplot as plt
from NEIVA.python_scripts.connect_with_mysql import*
from sqlalchemy import text

def fire_type():
    bk_db=connect_db('backend_db')
    dd=pd.read_sql(text('select * from bkdb_info_efcol'), con=bk_db)
    return list(dd['fire_type'].unique())

def table_info (database, fire_type):
    bk_db=connect_db('backend_db')
    dd=pd.read_sql(text('select * from bkdb_info_table_name'), con=bk_db)
    dd_final=dd[dd['db']==database].reset_index(drop=True)
    dd_final=dd_final[dd_final['fire_type'].str.contains(fire_type)]
    return dd_final[['tbl_name','measurement_type','fire_type','pollutant_category','study','source','doi']]
    
def summary_table (fire_type, measurement_type):
    bk_db=connect_db('backend_db')
    dd=pd.read_sql(text('select * from bkdb_info_efcol'), con=bk_db)
    dd=dd[dd['fire_type']==fire_type].reset_index(drop=True)
    if measurement_type=='all':
        dd=dd
    else:
        dd=dd[dd['measurement_type']==measurement_type].reset_index(drop=True)
        
    if fire_type == 'cookstove':
        return dd[['efcol','measurement_type','MCE','fuel_type','cookstove_name','cookstove_type','study']]
    else:
        return dd[['efcol','measurement_type','MCE','fuel_type','study']]

def display_pollutant_category():
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')
    
    rdf=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
    ll=rdf['pollutant_category'].unique()
    return list(ll)

def property_variables ():
    bk_db=connect_db('backend_db')
    df=pd.read_sql(text('select * from property_surrogate_info'), con=bk_db)
    return df[['column name', 'description', 'unit']]
    

def model_surrogates(chem):
   bk_db=connect_db('backend_db')
   output_db=connect_db('neiva_output_db')
   pp=pd.read_sql(text('select * from Property_Surrogate'), con=output_db)
   return  pp[chem][pp[chem].notna()].unique()

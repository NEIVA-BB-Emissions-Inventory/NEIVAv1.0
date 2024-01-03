#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  1 13:30:57 2024

@author: samiha
"""
import pandas as pd
import numpy as np
from sqlalchemy import text

from NEIVA.python_script.connect_with_mysql import *
legacy_db=connect_db('legacy_db')
raw_db=connect_db('raw_db')
primary_db=connect_db('primary_db')
bk_db=connect_db('backend_db')

def percent_c_correction_factor(dd, table_name):
    '''
    Apply carbon correction factor to emission factor data based on the specified table name.
    
    This function retrieves correction factors from the 'bkdb_correction_factor' database, 
    applies them to the appropriate columns of the provided dataset, and then returns the 
    corrected dataset. Specific corrections are applied to CO2 and CO based on the dataset's 
    origin.

    Parameters:
    - dd (DataFrame): Input dataset containing emission factor data.
    - table_name (str): Name of the table to determine which correction factors to use.

    Returns:
    - DataFrame: Corrected dataset with applied correction factors.
    '''
    
    cc=pd.read_sql(text('select * from bkdb_correction_factor'), con=bk_db)
    tbl_name = list(cc['pdb_table_name'][cc['pdb_table_name'].isin([table_name])].unique())
    if len(tbl_name)==1:
        columns =list(cc['efcol'][cc['pdb_table_name']==tbl_name[0]])
        cf = list(cc['correction factor'][cc['pdb_table_name']==tbl_name[0]])
        co2_ind = dd[dd['id']=='InChI=1S/CO2/c2-1-3'].index[0]
        co_ind = dd[dd['id']=='InChI=1S/CO/c1-2'].index[0]
        
        for i in range(len(columns)):
            dd[columns[i]]=dd[columns[i]]*cf[i]
            if columns[i].find('koss18')!=-1:
                dd.loc[co2_ind, columns[i]] = dd[columns[i]].iloc[co2_ind]/cf[i]
                dd.loc[co_ind, columns[i]] = dd[columns[i]].iloc[co_ind]/cf[i]
    return dd

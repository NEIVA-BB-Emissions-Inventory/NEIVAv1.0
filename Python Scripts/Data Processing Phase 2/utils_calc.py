#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 19:58:10 2022
@author: samiha
"""

import pandas as pd
import numpy as np

from utils import GrpCol

from connect_with_mysql import connect_db
bk_db=connect_db('backend_db')


def get_ind(intdf):
    com_df=pd.read_sql('select * from bkdb_imp_comList_for_calc', con=bk_db)
    co_ind=intdf[intdf['id']==com_df['id'][0]].index[0]
    co2_ind=intdf[intdf['id']==com_df['id'][2]].index[0]
    no_ind=intdf[intdf['id']==com_df['id'][1]].index[0]
    no2_ind=intdf[intdf['id']==com_df['id'][3]].index[0]
    return co_ind,co2_ind,no_ind,no2_ind

def calc_NOx_as_NO (intdf):
    '''
    Calculates the NOx values which is equivalent to NO in the input dataframe.
    This function fetches the values of NO and NO2 from the dataframe and 
    calculates their combined value as NOx equivalent to NO, using their 
    respective molecular weights for conversion.
    '''
    com_df=pd.read_sql('select * from bkdb_imp_comList_for_calc', con=bk_db)
    efcol=GrpCol(intdf)[2]
    
    no_ind=intdf[intdf['id']==com_df['id'].iloc[1]].index[0]
    no2_ind=intdf[intdf['id']==com_df['id'].iloc[3]].index[0]
    
    NOx_ind=intdf[intdf['id']=='NOx_as_NO'].index[0]
    
    for col in efcol:
        intdf.loc[NOx_ind, col] = intdf[col].iloc[no2_ind] * (30/46) + intdf[col].iloc[no_ind]
    
    return intdf





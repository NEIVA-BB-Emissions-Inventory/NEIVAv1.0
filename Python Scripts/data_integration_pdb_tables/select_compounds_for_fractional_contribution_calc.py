#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 12:47:02 2023

@author: samiha
"""
import pandas as pd
import numpy as np

from connect_with_mysql import *
legacy_db=connect_db('legacy_db')
raw_db=connect_db('raw_db')
primary_db=connect_db('primary_db')
bk_db=connect_db('backend_db')


def import_fc_dataset(nmogdf,lc_spec_df):
    '''
    Imports specific and simple fractional contribution datasets to the backend database.
    
    The function processes two dataframes, nmogdf and lc_spec_df, to derive two resulting datasets: 
    specific_fc_df and simple_fc. The specific dataset contains aligned lumped compounds and speciation,
    while the simple dataset consists of single lumped compounds for each formula.

    Parameters:
    - nmogdf: Dataframe with total NMOG data.
    - lc_spec_df: Dataframe with lumped compounds and their speciation.
    
    Returns:
    None. The function performs in-place modifications and exports results to the database.
    '''    
    # Construct specific fractional contribution dataset with aligned lumped compounds and speciation.
    specific_fc_df=lc_spec_df[GrpCol(lc_spec_df)[1]+['id','study']]
    specific_fc_df.to_sql(name='bkdb_fc_calc_specific', con=bk_db, if_exists='replace', index=False)
    
    # Extract simple fractional contribution dataset with only one lumped compound per formula.
    simple_fc=nmogdf[nmogdf['formula'].isin(GrpFormula(nmogdf)[3])]
    # Exclude entries present in specific_fc_df.
    simple_fc=simple_fc[~simple_fc['formula'].isin(specific_fc_df['formula'].tolist())]
    
    # Rearrange columns and export to backend database.
    simple_fc=simple_fc[GrpCol(simple_fc)[1]+['id','study']]
    simple_fc.to_sql(name='bkdb_fc_calc_simple',con=bk_db, if_exists='replace', index=False)
    return

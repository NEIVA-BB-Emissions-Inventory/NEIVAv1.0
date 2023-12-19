#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 19:40:04 2022
@author: Samiha Shahid
"""

import pandas as pd
import numpy as np
import pubchempy as pcp
from order_formula import exact_mm_calulator

'''
Establishing Database Connections:
This section imports the necessary functions to connect to various 
databases and then initializes connections to 
five specific databases: NEIVA_db, legacy_db, raw_db, primary_db, and backend_db.
'''
from connect_with_mysql import *
n_con=connect_db('NEIVA_db')
legacy_db=connect_db('legacy_db')
raw_db=connect_db('raw_db')
primary_db=connect_db('primary_db')
bk_db=connect_db('backend_db')
    
def AltName(df,df_altName):
    """
    Replaces the 'compound' names in the input dataframe 'df' based on 
    the 'altered_name' column of 'df_altName'. 

    Both 'df' and 'df_altName' are expected to have identical structures, 
    but 'df_altName' contains an additional 'altered_name' column.
    
    Parameters:
    - df: Dataframe containing original compound names.
    - df_altName: Dataframe with the same structure as 'df' but has an additional 
      'altered_name' column indicating the new names for compounds.

    Returns:
    - Updated dataframe with altered compound names.
    """
    iid=df_altName['id'][df_altName['altered_name'].notna()].to_list()
    altname=df_altName['altered_name'][df_altName['altered_name'].notna()].to_list()
    for i in range(len(iid)):
        iind=df[df['id']==iid[i]].index
        df.loc[iind,'compound']=altname[i]
    return df

def GrpCol(df):
    '''
    Categorizes columns of the dataframe based on their naming patterns and content.

    Parameters:
    - df: Dataframe to process.

    Returns:
    - arranged_col: Columns in the order of identity, emission factor, and ID.
    - identityCol: Columns identifying the chemical compounds.
    - efcol: Emission factor columns.
    - idcol: Column containing compound IDs.
    - avgcol: Columns with average values.
    - ncol: Columns indicating the number of data/studies.
    - stdcol: Columns showing standard deviations.
    '''
    identityCol=['mm','formula','compound','pollutant_category']
    efcol=df.filter(like='EF').columns.tolist()
    avgcol=df.filter(like='AVG').columns.tolist()
    ncol=df.filter(like='N').columns.tolist()
    stdcol=df.filter(like='STD').columns.tolist()
    idcol=['id']
    arranged_col=identityCol+efcol+idcol
    return arranged_col, identityCol, efcol, idcol, avgcol, ncol, stdcol


def sort_nmog(nmogdf):
    '''
    Sorts the NMOG (Non-Methane Organic Gases) dataframe based on molecular weight and 
    different criteria of the 'id' column.
    
    The function primarily sorts the dataframe first by the molecular mass. Within each molecular mass, 
    it further sorts rows based on the presence of 'InChI' in the 'id' column, followed by 
    rows without 'InChI' but not in the rdb_hatch15 dataset, and then rows with ids present 
    in the rdb_hatch15 dataset.
    
    Parameters:
    - nmogdf: Input dataframe containing NMOG data.
    
    Returns:
    - nmogdf_sorted: A sorted dataframe based on the above criteria.
    '''
    # Loadig the 'rdb_hatch15' dataste.
    hid = pd.read_sql('select * from rdb_hatch15',con=raw_db)
    hid = hid[hid['h_id'].notna()].reset_index(drop=True)
    
    # Sorting nmogdf by molecular mass.
    nmogdf=nmogdf.sort_values(by=['mm']).reset_index(drop=True)
    
    # List of unique formula of 'nmogdf'
    uf=list(nmogdf['formula'].unique())
    
    # Sorting each unique formula based on the criteria.
    nmogdf_sorted=pd.DataFrame()
    for f in uf:
        formula_df = nmogdf[nmogdf['formula']==f]
        # Select rows where 'id' contains 'InChI'
        aa1=formula_df[formula_df['id'].str.contains('InChI')]
        aa1=aa1.sort_values(by='id')
        
        # Select rows where 'id' does not contains InChI and is not in hid.
        aa2=formula_df[~formula_df['id'].str.contains('InChI')][~formula_df['id'].isin(hid['h_id'])]
        
        # Select rows where 'id' does not contains InChI but is in hid.
        aa3=formula_df[~formula_df['id'].str.contains('InChI')][formula_df['id'].isin(hid['h_id'])]
        
        aa=aa2.append(aa1).append(aa3)
        nmogdf_sorted=nmogdf_sorted.append(aa).reset_index(drop=True)
    
    return nmogdf_sorted

def rearrange_col_finaldf (df):
    '''
    Rearranges the columns of the input dataframe based on a specific order 
    derived from an external 'bkdb_info_efcol' in the database.
    
    The function fetches fire types from the 'bkdb_info_efcol', and then for each 
    fire type, it looks for matching columns in the input dataframe. Finally, 
    columns are rearranged based on this order, along with a few specific columns 
    from the input dataframe.
    
    Parameters:
    - df: Input dataframe whose columns need to be rearranged.
    
    Returns:
    - df: A dataframe with columns rearranged in the desired order.
    '''
    
    efcoldf=pd.read_sql('select * from bkdb_info_efcol', con=bk_db)
    ft_ll=efcoldf['fire_type'].unique()
    
    col_ll=[]
    for i in range(len(ft_ll)):
        col=(list(df.columns[df.columns.str.contains(ft_ll[i].replace(' ','_'))]))
        for k in range(len(col)):
            col_ll.append(col[k])
    
    total_col=GrpCol(df)[1]+col_ll+['id']
    df=df[total_col]
    return df


def str_float(df, col):
    '''
    Converts valid string representations of numbers in a specified column of a dataframe 
    into floating-point numbers. The function ignores NaN and None values.
    
    Parameters:
        df (pd.DataFrame): The input dataframe that contains the column to be processed.
        col (str): The column name which contains string representations of numbers 
                   that need to be converted to float.
    
    Returns:
        pd.DataFrame: Updated dataframe with the specified column's string values 
                      converted to float.
    '''
    for i in range(len(df)):
        if str(df[col].iloc[i])!='nan':
            if str(df[col].iloc[i])!='None':
                df.loc[i,col]=float(df[col].iloc[i])
    return df





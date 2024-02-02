#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 19:40:04 2022
@author: Samiha Shahid
"""

import pandas as pd
import numpy as np
from sqlalchemy import text
from NEIVA.python_scripts.data_integration_process.sort_molec_formula import exact_mm_calulator

'''
Establishing Database Connections:
This section imports the necessary functions to connect to various 
databases and then initializes connections to 
five specific databases: NEIVA_db, legacy_db, raw_db, primary_db, and backend_db.
'''
from NEIVA.python_scripts.connect_with_mysql import connect_db
# n_con=connect_db('NEIVA_db')
# legacy_db=connect_db('legacy_db')
# raw_db=connect_db('raw_db')
# primary_db=connect_db('primary_db')
# bk_db=connect_db('backend_db')

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
    
    bk_db=connect_db('backend_db')
    efcoldf=pd.read_sql(text('select * from bkdb_info_efcol'), con=bk_db)
    ft_ll=efcoldf['fire_type'].unique()
    
    col_ll=[]
    for i in range(len(ft_ll)):
        col=(list(df.columns[df.columns.str.contains(ft_ll[i].replace(' ','_'))]))
        for k in range(len(col)):
            col_ll.append(col[k])
    
    total_col=GrpCol(df)[1]+col_ll+['id']
    df=df[total_col]
    return df

def assign_study_column(nmogdf):
    '''
    Assigns appropriate study names to each row in the 'nmogdf' dataframe based on 
    the non-null emission factor columns present in that row.
    
    The function fetches the study info from 'bkdb_info_efcol' table, reorders 
    them based on year, and then maps the study name to the corresponding row in 
    'nmogdf' dataframe. In case of overlapping studies, specific naming adjustments 
    are made.
    
    Parameters:
        nmogdf (pd.DataFrame): The dataframe containing emission factor data.
    
    Returns:
        pd.DataFrame: The dataframe 'nmogdf' updated with 'study' column.
    '''
    bk_db=connect_db('backend_db')
    efcoldf=pd.read_sql(text('select * from bkdb_info_efcol'), con=bk_db)
    efcols=GrpCol(nmogdf)[2]
    
    efcoldf=efcoldf.sort_values(by=['year','year_akagi_data'], ascending=False)
    efcoldf=efcoldf.reset_index(drop=True)
    
    # overlapping study between ldb and rdb
    for i in efcoldf[efcoldf['study']=='Akagi_11(stockwell15)'].index:
        efcoldf.loc[i,'study']='stockwell15'
    
    for i in range(len(nmogdf)):
        efdat=list(nmogdf[efcols].iloc[i].dropna().index)
        st=list(efcoldf['study'][efcoldf['efcol'].isin(efdat)].unique())
        nmogdf.loc[i,'study']=(',').join(st)
    return nmogdf

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

def merge_pm(mdf):
    ind=list(mdf[mdf['pollutant_category']=='PM total']\
        [mdf['compound'].str.contains('PM')]\
        [~mdf['id'].isin(['PM10','PM2.5_npb','PM2.5_pb','PM2.5_ipcc'])].index)
        
    efcol=list(mdf.columns[mdf.columns.str.contains('EF')])
    vals=list(mdf[efcol][mdf.index.isin(ind)].mean().values)
    pmind=mdf[mdf['id']=='PM2.5'].index[0]
    for i in range(len(efcol)):
        mdf.loc[pmind,efcol[i]]=vals[i]
    
    dropind=set(ind)-{pmind}
    mdf=mdf.drop(index=dropind)
    mdf=mdf.reset_index(drop=True)
    return mdf


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

def GrpFormula(df):
    '''
    Categorizes chemical formulas from the input dataframe.
    
    The categorization is done based on formula occurrence and the presence of lumped compounds.
    A lumped compound is determined if its ID doesn't match the pdb_hatch15 table or an InChI structure.

    Parameters:
    - df: Dataframe containing chemical data.

    Returns:
    - f_no_spec: Formulas without speciation.
    - f_spec_no_lc: Formulas with speciation but without lumped compounds.
    - f_spec_multiple_lc: Formulas with speciation and multiple lumped compounds.
    - f_spec_one_lc: Formulas with speciation and a single lumped compound.
    
    '''
    # Loading the hatch15 from primary database.
    hid = pd.read_sql('select * from pdb_hatch15',con=primary_db)
    
    # Unique formula of 'df' dataframe.
    all_formula=df['formula'].unique().tolist()
    
    # Identify formula that only appear once in the dataframe.
    f_no_spec=[]
    for f in all_formula:
        if len(df[df['formula']==f])==1: 
            f_no_spec.append(f)
    
    # Formula with more than one occurance.
    f_spec=set(all_formula)-set(f_no_spec)
    
    # Formulas with speciation and without lumped compounds.
    f_spec_no_lc=[]
    for f in f_spec:
        aa=df[df['formula']==f]
        bb=aa[aa['id'].str.contains('InChI',na=False)]
        if len(aa)==len(bb):
            f_spec_no_lc.append(f)
    
    # Formulas that are not in 'f_spec_no_lc'
    f_spec_lc=f_spec-set(f_spec_no_lc)
    
    # Formulas with multiple lumped compound.
    f_spec_multiple_lc=[]
    for f in f_spec_lc:
        aa=df[df['formula']==f]
        aa=aa[~aa['id'].str.contains('InChI',na=False)]
        aa=aa[~aa['id'].isin(hid['id'].tolist())]
        if len(aa)>1:
             f_spec_multiple_lc.append(f)
    # Formulas with a single lumped compound.
    f_spec_one_lc=[]
    for f in f_spec_lc:
        aa=df[df['formula']==f]
        aa=aa[~aa['id'].str.contains('InChI',na=False)]
        aa=aa[~aa['id'].isin(hid['id'])]
        if len(aa)==1:
            f_spec_one_lc.append(f)
    
    return f_no_spec, f_spec, f_spec_multiple_lc, f_spec_one_lc 

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
    
    efcoldf=pd.read_sql('select * from bkdb_info_efcol', con=bk_db)
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
    
    cc=pd.read_sql('select * from bkdb_correction_factor', con=bk_db)
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




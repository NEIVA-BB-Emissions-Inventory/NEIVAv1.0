#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 12:38:36 2022
@author: Samiha Shahid
"""
import pandas as pd
from sqlalchemy import text
from sqlalchemy import create_engine
# import numpy as np
# import pubchempy as pcp
from NEIVA.python_scripts.data_integration_process.data_formatting_functions import AltName,GrpCol
from NEIVA.python_scripts.data_integration_process.categorize_chemical_formula import *

'''
Establishing Database Connections:
This section imports the necessary functions to connect to various 
databases and then initializes connections to 
five specific databases: NEIVA_db, legacy_db, raw_db, primary_db, and backend_db.
'''
from NEIVA.python_scripts.connect_with_mysql import connect_db
# legacy_db=connect_db('legacy_db')
# raw_db=connect_db('raw_db')
# primary_db=connect_db('primary_db')
# bk_db=connect_db('backend_db')


def get_lumped_com_id_df(f_spec_lc,df):
    '''
    Extract lumped compound IDs for a given formula list from a dataframe.
    
    Parameters:
    - f_spec_lc: List of formulas to filter.
    - df: Dataframe to extract IDs from.
    
    Process:
    1. Exclude hatch15 IDs.
    2. Filter for rows based on the provided formula list.
    3. Remove rows with IDs containing 'InChI'.
    4. Exclude rows that match hatch15 isomers.
    
    Returns:
    - Dataframe containing lumped compound IDs that don't have an InChI.
    '''
    raw_db=connect_db('raw_db')
    hid = pd.read_sql(text('select * from rdb_hatch15'),con=raw_db)
    hid = hid[hid['h_id'].notna()].reset_index(drop=True) # exclude h15 isomers
    
    df=df[df['formula'].isin(f_spec_lc)].reset_index(drop=True)
    df=df[~df['id'].str.contains('InChI')]
    df=df[~df['id'].isin(hid['h_id'])]
    df=df.reset_index(drop=True)
    return df
        
def select_max_len_CompoundStr(com_ll,id_ll):
    '''
    Given a list of compounds and their IDs, return the ID of the compound 
    with the longest name.

    Parameters:
    - com_ll: List of compound names.
    - id_ll: Corresponding list of compound IDs.

    Returns:
    - ID of the compound with the longest name.
    '''
    
    len_ll=[]
    for com in com_ll:
        len_ll.append(len(com))
    selected_ind=len_ll.index(max(len_ll))
    selected_id=id_ll[selected_ind]
    return  selected_id  

def select_id_df(iddf,formula):
    '''
    For each formula in the provided list, the function selects the compound
    with the longest name among those with multiple lumped compounds. 
    It then returns a dataframe containing only the selected compounds.

    Parameters:
    - iddf: Dataframe containing compound details.
    - formula: List of formulae to consider.

    Returns:
    - A dataframe with compounds having the longest names for each formula.
    '''    
    slc_id=[]
    for f in formula:
        com_ll=iddf['compound'][iddf['formula']==f].tolist()
        id_ll=iddf['id'][iddf['formula']==f].tolist()
        sl_id=select_max_len_CompoundStr(com_ll,id_ll)
        slc_id.append(sl_id)
    iddf=iddf[iddf['id'].isin(slc_id)]
    iddf=iddf.reset_index(drop=True)
    return iddf
    
def merge_rows(df):
    '''
    Transposes multiple rows with the same formula into a single row 
    by taking the mean of their values.

    Parameters:
    - df: Dataframe containing multiple rows for specific formulas.

    Returns:
    - A dataframe with each unique formula represented by a single row.
    '''    
    uf=df['formula'].unique().tolist()
    efcols=GrpCol(df)[2]
    
    ef_df=pd.DataFrame()
    for formula in uf:
        aa=pd.DataFrame(df[efcols][df['formula']==formula].mean()).transpose()
        # ef_df=ef_df.append(aa)
        ef_df = pd.concat([ef_df,aa], ignore_index=True)
    ef_df=ef_df.reset_index(drop=True)
    ef_df['formula']=uf
    return ef_df

def check_r_iddf (iddf, r_iddf):
    '''
    Validates the results of the merge_rows() function by comparing the number of columns 
    in the original and reduced dataframes for each unique formula.

    Parameters:
    - iddf: Original dataframe containing multiple rows for specific formulas.
    - r_iddf: Reduced dataframe after applying the merge_rows() function.

    Returns:
    - "Clear" if the validation passes for all formulas, otherwise prints the formulas that didn't match.
    '''    
    uf=list(iddf['formula'].unique())
    ch=[]
    
    for u in uf:
        ddf=ddf=iddf[iddf['formula']==u].reset_index(drop=True)
        
        col=[]
        for i in range(len(ddf)):
            xx=set(ddf.iloc[i].dropna().index)-{'mm','formula','compound','pollutant_category','id','study'}
            col=col+list(xx)
            
        aa=len(set(col))
        bb=len(set(r_iddf[r_iddf['formula']==u].T.dropna().index)-{'mm','formula','compound','pollutant_category','id','study'})
            
        if aa==bb:
            print('______________________________________________________________________________________________________________')
            print('Formula-'+u, 'Merged compounds-')
            print(iddf[['compound','study']][iddf['formula']==u].reset_index(drop=True))
            print('______________________________________________________________________________________________________________')
        else:
            print(u, 'error: Check the merging step')
            ch.append(u)
    if len(ch)==0:
        print('Rows merging process is verified.')
        return 'Clear'
                                
def alter_name_slc_iddf (slc_iddf):
    
    com_name = ['Methyl Cyclopentadiene ( isomer 1, C6H8)',
                 'Other C6H10 (isomer_1)',
                 'Hexenes (sum of 3 isomers)',
                 'assorted amides',
                 '2+3-methylpentane',
                 'assorted amines',
                 'assorted hcs',
                 'undaturated C6 cyclic carboxylic acid',
                 'terpenes(-pinene)']
    
    altered_name = ['Methyl Cyclopentadiene isomers',
                     'C6H10 isomers',
                     'Hexenes+C6H12 isomers',
                     'C4H9NO amides',
                     '2-methylpentane+3-methylpentane',
                     'C4H11NO amines',
                     'C7H12 isomers',
                     'unsaturated C6 cyclic carboxylic acid',
                     'Monoterpenes']
    
    for i in range(len(com_name)):
        ind=slc_iddf[slc_iddf['compound']==com_name[i]].index[0]
        slc_iddf.loc[ind,'compound']=altered_name[i]
    
    return slc_iddf


def merge_lumped_compound_same_formula(nmogdf):
    '''
    Merges rows with the same chemical formula if they represent lumped compounds.
    A limped compound is identified when its ID is not derived from retention index 
    (pdb_hatch15) and an InChI.
    
    Input:
        - NMOC_g dataset
    
    Returns:
        - A consolidated dataframe of lumped compounds (r_iddf)
        - IDs of lumped compounds (iddf)    
    '''
    
    # List of chemical formula that has multiple lumped compounds.
    f_spec_multiple_lc=assign_formula_type(nmogdf)[2]
    
    # Dataframe with the lumped compound ids
    iddf=get_lumped_com_id_df(f_spec_multiple_lc,nmogdf)
    
    #Improting to backend_db
    engine = create_engine("mysql+pymysql://root:root@localhost/"+'backend_db')
    iddf[GrpCol(iddf)[1]+['id']].to_sql(name='bkdb_nmog_MultLumCom',con=engine, if_exists='replace', index=False)
    
    # List of formula where the multiple compounds remain unmerged
    slcdf=['C4H6'] 
    
    # Excluding the formulas where lumped compounds are not merged.
    f_spec_multiple_lc=set(f_spec_multiple_lc)-set(slcdf) 
    
    # Updating the 'iddf' dataframe with the new formula list.
    iddf=get_lumped_com_id_df(f_spec_multiple_lc,nmogdf) 
    
    # Selects the identifier of the lumped compound for a specific formula.
    # The compound name with the longest string length is chosen. 
    slc_iddf=select_id_df(iddf,f_spec_multiple_lc) 
    
    # Store the identity columns from 'slc_iddf'.
    slc_iddf=slc_iddf[GrpCol(slc_iddf)[1]+['id']] 
    
    # Save the 'slc_iddf' to backend_db for user review. 
    slc_iddf.to_sql(name='bkdb_nmog_MultLumCom_slc_id',con=engine, if_exists='replace',index=False) 
    
    # Changing compound names
    slc_iddf=alter_name_slc_iddf (slc_iddf)
    
    # Transpose the EF columns so that multiple rows of a single formula become a single row.
    iddf_ef=merge_rows(iddf) 
    
    r_iddf=pd.DataFrame()
    r_iddf=slc_iddf.merge(iddf_ef, on='formula', how='left')
    
    r_iddf=r_iddf[GrpCol(r_iddf)[0]] # the columns are arranged
    
    check_r_iddf (iddf, r_iddf)
    print('+--------------------------------------------------------------------------------------------+')
    print('Length of merged dataframe : '+str(len(r_iddf)))
    print('Length of dataframe containing lumped compound(>1) of a unique formula : '+str(len(iddf)))
    print('+--------------------------------------------------------------------------------------------+')
    return r_iddf, iddf    

def insert_rdf_nmogdf(nmogdf,rdf,df):
    '''
    Function to Update NMOC_g Dataset:
    This function appends the reduced dataframe 'rdf' to 'nmogdf' after removing rows from 'nmogdf' that are already present in 'df'.
    After the insertion, it prints the updated length of the NMOC_g dataset.
    Input: 
    - nmogdf: Original NMOC_g dataset
    - rdf: Reduced dataframe to be appended
    - df: Reference dataframe
    Returns:
    - Updated NMOC_g dataset
    '''
    nmogdf=nmogdf[~nmogdf['id'].isin(df['id'])]
    nmogdf = pd.concat([nmogdf,rdf], ignore_index=True)
    # nmogdf=nmogdf.append(rdf).reset_index(drop=True)
    print('Length of NMOC_g dataset: '+str(len(nmogdf)))
    return nmogdf


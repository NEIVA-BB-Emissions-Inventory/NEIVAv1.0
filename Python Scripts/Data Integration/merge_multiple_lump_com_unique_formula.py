#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 12:38:36 2022
@author: Samiha Shahid
"""
import pandas as pd
import numpy as np
import pubchempy as pcp
from utils import GrpFormula,AltName,GrpCol

'''
Establishing Database Connections:
This section imports the necessary functions to connect to various 
databases and then initializes connections to 
five specific databases: NEIVA_db, legacy_db, raw_db, primary_db, and backend_db.
'''
from connect_with_mysql import*
n_con=connect_db('NEIVA_db')
legacy_db=connect_db('legacy_db')
raw_db=connect_db('raw_db')
primary_db=connect_db('primary_db')
bk_db=connect_db('backend_db')


def Get_lcid_df(f_spec_lc,df):
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
    
    hid = pd.read_sql('select * from rdb_hatch15',con=raw_db)
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

def slc_iddf_Mult_lcid(iddf,formula):
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
    
def Mult_row2single_Row(df):
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
        ef_df=ef_df.append(aa)
    ef_df=ef_df.reset_index(drop=True)
    ef_df['formula']=uf
    return ef_df

def check_r_iddf (iddf, r_iddf):
    '''
    Validates the results of the Mult_row2single_Row() function by comparing the number of columns 
    in the original and reduced dataframes for each unique formula.

    Parameters:
    - iddf: Original dataframe containing multiple rows for specific formulas.
    - r_iddf: Reduced dataframe after applying the Mult_row2single_Row() function.

    Returns:
    - "Clear" if the validation passes for all formulas, otherwise prints the formulas that didn't match.
    '''    
    uf=list(iddf['formula'].unique())
    ch=[]
    
    for u in uf:
        ddf=ddf=iddf[iddf['formula']==u].reset_index(drop=True)
        
        col=[]
        for i in range(len(ddf)):
            xx=set(ddf.iloc[i].dropna().index)-{'mm','formula','compound','pollutant_category','id'}
            col=col+list(xx)
            
        aa=len(set(col))
        bb=len(set(r_iddf[r_iddf['formula']==u].T.dropna().index)-{'mm','formula','compound','pollutant_category','id'})
            
        if aa==bb:
            print('y')
        else:
            print(u, 'n')
            ch.append(u)
    if len(ch)==0:
        return 'Clear'
                        


def reduce_multiple_LumCom(nmogdf):
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
    f_spec_multiple_lc=GrpFormula(nmogdf)[2]
    
    # Dataframe with the lumped compound ids
    iddf=Get_lcid_df(f_spec_multiple_lc,nmogdf)
    
    #Improting to backend_db
    iddf[GrpCol(iddf)[1]].to_sql(name='bkdb_nmog_MultLumCom',con=bk_db, if_exists='replace', index=False)
    
    # List of formula where the multiple compounds remain unmerged
    slcdf=['C4H6'] 
    
    # Excluding the formulas where lumped compounds are not merged.
    f_spec_multiple_lc=set(f_spec_multiple_lc)-set(slcdf) 
    
    # Updating the 'iddf' dataframe with the new formula list.
    iddf=Get_lcid_df(f_spec_multiple_lc,nmogdf) 
    
    # Selects the identifier of the lumped compound for a specific formula.
    # The compound name with the longest string length is chosen. 
    slc_iddf=slc_iddf_Mult_lcid(iddf,f_spec_multiple_lc) 
    
    # Store the identity columns from 'slc_iddf'.
    slc_iddf=slc_iddf[GrpCol(slc_iddf)[1]+['id']] 
    
    # Save the 'slc_iddf' to backend_db for user review. 
    slc_iddf.to_sql(name='bkdb_nmog_MultLumCom_slc_id',con=bk_db, if_exists='replace',index=False) 
    
    # Load the 'bkdb_nmog_MultLumCom_slc_id_altName' which is similar to 'bkdb_nmog_MultLumCom_slc_id' 
    # but with an added 'altered_name' column. This column has modified compound names from the 'compound' column.
    df_altName = pd.read_sql('select * from bkdb_nmog_MultLumCom_slc_id_altName', con=bk_db)
    slc_iddf = AltName(slc_iddf, df_altName)
    
    # Transpose the EF columns so that multiple rows of a single formula become a single row.
    iddf_ef=Mult_row2single_Row(iddf) 
    
    r_iddf=pd.DataFrame()
    r_iddf=slc_iddf.merge(iddf_ef, on='formula', how='left')
    
    r_iddf=r_iddf[GrpCol(r_iddf)[0]] # the columns are arranged
    
    check_r_iddf (iddf, r_iddf)
    print('Length of r_iddf: '+str(len(r_iddf)))
    print('Length of iddf: '+str(len(iddf)))
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
    nmogdf=nmogdf.append(rdf).reset_index(drop=True)
    print('Length of NMOC_g dataset: '+str(len(nmogdf)))
    return nmogdf


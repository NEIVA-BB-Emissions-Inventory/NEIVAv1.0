#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 20:02:38 2022

@author: samiha
"""
import pandas as pd
import numpy as np
from sqlalchemy import text

from NEIVA.python_scripts.data_integration_process.data_formatting_functions import GrpCol
from NEIVA.python_scripts.data_integration_process.categorize_chemical_formula import assign_formula_type
from NEIVA.python_scripts.data_processing_steps.data_calculations import get_ind

from NEIVA.python_scripts.connect_with_mysql import *



def test(lumpval,spec_ef_sum):
    '''
    Tests if 'spec_ef_sum' falls within a range around 'lumpval'.
    Parameters:
    - lumpval: The reference value.
    - spec_ef_sum: The value to be tested.
    Returns:
    - 'yes' if 'spec_ef_sum' is within 2.5 times 'lumpval' of the range centered around 'lumpval'.
    - 'no' otherwise.
    '''
    # Define maximum and minimum values based on 'lumpval'
    max_val=lumpval+2.5*lumpval
    min_val=lumpval-2.5*lumpval
   
    # Ensure that 'min_val' is not negative.
    if min_val<0:
        min_val=0
    
    # Check id 'spec_ef_sum' is within the specified range.     
    if (spec_ef_sum<=max_val) and (spec_ef_sum>min_val) :
        return 'yes'
    else:
        return 'no'

def get_ind_list_sfc(sfc,formula):
    '''
    Returns indices of lumped ids within a formula of a DataFrame.
    '''
    
    # Select rows with the specified formula from the 'sfc' dataframe.
    df=sfc[sfc['formula']==formula]
    # Identify indices of lumped compounds within the DataFrame.
    lump_inds=df[~df['id'].str.contains('InChI')].index.tolist()
    
    # Add the length of the DataFrame.
    lump_inds.append(len(df))
    
    ind_ll=[]
    for i in range(len(lump_inds)-1):
        ll=(np.arange(lump_inds[i],lump_inds[i+1])).tolist()
        ind_ll.append(ll)
    return ind_ll

def Get_fc_calc(fc,df):
    '''
    Calculates fractional contributions for average columns in the 'df' DataFrame based on 'fc' DataFrame.

    Parameters:
    - fc: The DataFrame containing fractional contribution data.
    - df: The DataFrame with emission data and average columns.

    Returns:
    - The modified 'df' DataFrame with fractional contribution adjustments.
    '''
    primary_db=connect_db('primary_db')

    # Load the hatch15 dataset.
    hid=pd.read_sql(text('select * from pdb_hatch15'), con=primary_db)    
    # Get a list of unique formulas from the 'fc' DataFrame.
    uf=fc['formula'].unique().tolist()
    # Get the list of average columns in 'df'.
    avgcols=GrpCol(df)[4]
    
    # Iterating over unique formulas
    for i in range(len(uf)):
        formula=uf[i]
        # Selcet the lumped ID of the specific formula.
        lumpid=fc['id'][fc['formula']==uf[i]][~fc['id'].str.contains('InChI')][~fc['id'].isin(hid['id'])].values[0]
        # Select speciation IDs for the specific formula.
        specids=fc['id'][fc['formula']==uf[i]][~fc['id'].isin([lumpid])].tolist()
        
        # Iterate over average columns.
        for col in avgcols:
            # get the value of the lumped ID.
            lumpval=df[col][df['id']==lumpid].values[0]
            # Sum the EF values of speciation compounds.
            spec_ef_sum=df[col][df['id'].isin(specids)].sum() # sum the ef of speciation.
            
            # Get the index of the lumped compounds.
            lumpind=df[df['id']==lumpid].index[0]
            # Get the index of the speciation compound.
            spec_inds=df[df['id'].isin(specids)].index.tolist()
            
            # Test if the sum of speciation EF values is within 2.5 times the lumped value.
            if test(lumpval,spec_ef_sum)=='yes':
                # Alpply fractional contribution to speciation indices.
                for ind in spec_inds: 
                    fc_val=df[col].iloc[ind]/spec_ef_sum
                    lump_ef_contribution=df[col].iloc[lumpind]*fc_val
                    ncol='N_'+col.split('AVG_')[1]
                    n=df[ncol].iloc[ind]
                    df.loc[ind,col]=(lump_ef_contribution+df[col].iloc[ind]*n)/(n+1)
                    
                    # if str(df[col].iloc[ind])!='nan':
                    #     df.loc[ind,ncol]=n+1
                
                df.loc[lumpind,col]=np.nan
                df.loc[lumpind,ncol]=np.nan
                df.loc[lumpind,'STD_'+col.split('AVG_')[1]]=0
    return df

def assign_fractional_contribution (df):
    '''
    Calculates fractional contribution (FC) values based on two fractional contribution datasets.
    Updates the 'df' DataFrame with FC values.
    
    Parameters:
    - df: The DataFrame containing emission data.
    Returns:
    - The 'df' DataFrame with FC values.
    '''
    bk_db=connect_db('backend_db')

    # Loading two FC datasets
    fc=pd.read_sql(text('select * from bkdb_fc_calc_simple'), con=bk_db) # Simple fractional contribution dataset
    sfc=pd.read_sql(text('select * from bkdb_fc_calc_specific'), con=bk_db) # Specific fractional contribution dataset
    
    # Separate lumped compounds from specific FC dataset 
    f_m_lids=assign_formula_type(sfc)[2] # Get formulas with multiple lumped IDs from 'sfc'
    fc2=sfc[~sfc['formula'].isin(f_m_lids)] # Dataframe without multiple lumped compounds.
    sfc=sfc[sfc['formula'].isin(f_m_lids)].reset_index(drop=True) # Assign multiple lumped compounds to 'sfc'
    
    # Append the seperated data back to the simple FC dataset.
    fc = pd.concat([fc, fc2], ignore_index=True)    
    # Get average columns in the in df.
    avgcols=GrpCol(df)[4] 
    
    # Calculate fractional contribution and update 'df'.
    df=Get_fc_calc(fc,df)
    
    # Get unique formula from 'sfc' dataframe.
    sfc_uf=sfc['formula'].unique().tolist() 
    
    # Iterate over unique formula of 'sfc'.
    for i in range(len(sfc_uf)):
        # Get indices of rows with the current formula.
        formula=sfc_uf[i]
        ll=get_ind_list_sfc(sfc,formula)
        # Iterate over the indices and calculate FC.
        for k in range(len(ll)):
            sdf=sfc[sfc.index.isin(ll[k])]
            df=Get_fc_calc(sdf,df)
    
    # Remove rows with misisng values of average EF columns.
    d_ind=set(df.index)-set(df[GrpCol(df)[4]].dropna(how='all').index)
    df=df.drop(index=d_ind)
    df=df.reset_index(drop=True)
    return df

def round_avg_cols (df):
     avgcols=GrpCol(df)[4] # Get the average EF columns.
     stdcols=GrpCol(df)[6] # Get the STD columns.
     
     for col in avgcols:
         df[col]=round(df[col],4)
     for col in stdcols:
         df[col]=round(df[col],4)
     return df



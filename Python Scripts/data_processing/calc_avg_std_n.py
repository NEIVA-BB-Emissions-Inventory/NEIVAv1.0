#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 12:17:34 2023

@author: samiha
"""

import pandas as pd
import numpy as np

from utils import GrpCol, GrpFormula
from utils_calc import get_ind

from connect_with_mysql import*
bk_db=connect_db('backend_db')
primary_db=connect_db('primary_db')

def assign_n_cols(df,efcoldf):
    '''
     Adds a 'N_fire_type' column to the 'df' DataFrame, showing the count of non-null values
    in certain columns for each unique 'fire_type' found in the 'efcoldf' DataFrame.
    '''
    # Get a list of unique 'fire type' values from 'efcoldf'
    ft_ll=efcoldf['fire_type'].unique().tolist() 
    # Loop through each unique fire_type
    for fire_type in  ft_ll:
        # Get a list of EF columns specific to the current 'fire type'.
        efcols= efcoldf['efcol'][efcoldf['fire_type']==fire_type].tolist()
        # Create a DataFrame containing only the EF columns of the current 'fire_type'.
        ft_df=df[efcols]
        # Assign a name for the 'N' column based on the fire_type.
        ncolName='N_'+fire_type.replace(' ','_')
        # Calculate and assign the count of non-null values to the 'N' column for each row.
        for i in range(len(ft_df)):
            df.loc[i,ncolName]=len(df[efcols].iloc[i].dropna().values)
    return df


def get_avg_df(df,efcoldf):
    '''
    Calculates average and standard deviation values for specific columns in the 'df' DataFrame, 
    grouped by unique 'fire_type' categories from 'efcoldf'. Returns a new DataFrame with the computed values.
    
    Parameters:
    - df: The DataFrame containing emission data.
    - efcoldf: The DataFrame with emission factor (EF) information and 'fire_type' categories.

    Returns:
    - A new DataFrame with average and standard deviation values for each 'fire_type' category.
    
    '''
    # Getting a list of unique fire type values.
    ft_ll=efcoldf['fire_type'].unique().tolist() # getting the fire type unique list.
    # Loop through each unique 'fire type'.
    for fire_type in  ft_ll:
        # Get a list of EF columns specific tot he current 'fire_type'.
        efcols= efcoldf['efcol'][efcoldf['fire_type']==fire_type].tolist()
        
        # Assign column names for average and standard deviation.
        colName='AVG_'+fire_type.replace(' ','_')
        stdcolName='STD_'+fire_type.replace(' ','_')
        
        # Calculate the average EF and assign to the new column.
        df[colName]=df[efcols].T.mean()
        # Calculate the standard deviation of the EF and assign to the new column.
        df[stdcolName]=df[efcols].T.std() 
    
    # Rearange columns in the new dataframe.
    col_ll=GrpCol(df)[1]+GrpCol(df)[4]+GrpCol(df)[5]+GrpCol(df)[6]+['id'] 
    avgdf=df[col_ll]
    
    # Calculate PM<2.5 values and replace PM2.5 with PM<2.5.
    
    avgcol=GrpCol(avgdf)[4]
    pmvals=avgdf[avgcol][avgdf['pollutant_category']=='PM total'][avgdf['compound'].str.contains('PM',na=False)].mean().values.tolist()
    pmind=avgdf[avgdf['id']=='PM2.5'].index[0]
    avgdf.loc[pmind,'compound']='PM<2.5'
    avgdf.loc[pmind,'id']='PM<2.5'
    
    # Assigning PM<2.5 values tp all average columns.
    for i in range(len(avgcol)):
        avgdf.loc[pmind,avgcol[i]]=pmvals[i]
    
    # Remove rows with total PM and reset the index.    
    total_pm_ind=set(avgdf[avgdf['pollutant_category']=='PM total'][avgdf['compound'].str.contains('PM',na=False)].index)-{pmind}
    
    avgdf=avgdf.drop(index=total_pm_ind)
    avgdf=avgdf.reset_index(drop=True)
    
    return avgdf

def round_avg_cols (df):
     avgcols=GrpCol(df)[4] # Get the average EF columns.
     stdcols=GrpCol(df)[6] # Get the STD columns.
     
     for col in avgcols:
         df[col]=round(df[col],4)
     for col in stdcols:
         df[col]=round(df[col],4)
     return df

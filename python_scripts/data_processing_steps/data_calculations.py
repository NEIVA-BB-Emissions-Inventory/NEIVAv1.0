#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 19:58:10 2022
@author: samiha
"""

import pandas as pd
import numpy as np
from sqlalchemy import text


from NEIVA.python_scripts.data_integration_process.data_formatting_functions import GrpCol

from NEIVA.python_scripts.connect_with_mysql import *


def get_ind(intdf):
    bk_db=connect_db('backend_db')

    com_df=pd.read_sql(text('select * from bkdb_imp_comList_for_calc'), con=bk_db)
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
    bk_db=connect_db('backend_db')

    com_df=pd.read_sql(text('select * from bkdb_imp_comList_for_calc'), con=bk_db)
    efcol=GrpCol(intdf)[2]
    
    no_ind=intdf[intdf['id']==com_df['id'].iloc[1]].index[0]
    no2_ind=intdf[intdf['id']==com_df['id'].iloc[3]].index[0]
    
    NOx_ind=intdf[intdf['id']=='NOx_as_NO'].index[0]
    
    for col in efcol:
        if str(intdf[col].iloc[NOx_ind])=='nan':
            intdf.loc[NOx_ind, col] = intdf[col].iloc[no2_ind] * (30/46) + intdf[col].iloc[no_ind]
    return intdf


def assign_data_count_column(df,efcoldf):
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


def calculate_average_fire_types(df,efcoldf):
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
    pmvals=avgdf[avgcol][avgdf['pollutant_category']=='PM total'][avgdf['compound'].str.contains('PM',na=False)][avgdf['id']!='PM10'][avgdf['id']!='PM2.5_ipcc'].mean().values.tolist()
    pmind=avgdf[avgdf['id']=='PM2.5'].index[0]
    avgdf.loc[pmind,'compound']='PM2.5*'
    avgdf.loc[pmind,'id']='PM2.5*'
    
    # Assigning PM<2.5 values tp all average columns.
    for i in range(len(avgcol)):
        avgdf.loc[pmind,avgcol[i]]=pmvals[i]
    
    # Remove rows with total PM and reset the index.    
    total_pm_ind=set(avgdf[avgdf['pollutant_category']=='PM total'][avgdf['compound'].str.contains('PM',na=False)][avgdf['id']!='PM10'].index)-{pmind}
    
    avgdf=avgdf.drop(index=total_pm_ind)
    avgdf=avgdf.reset_index(drop=True)
    
    return avgdf




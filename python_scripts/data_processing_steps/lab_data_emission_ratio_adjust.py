#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 19:58:18 2022

@author: samiha
"""
import pandas as pd
import numpy as np
from sqlalchemy import text

from NEIVA.python_scripts.data_processing_steps.data_calculations import get_ind
from NEIVA.python_scripts.data_processing_steps.info_table_sorting_functions import sort_by_fire_type_year
from NEIVA.python_scripts.connect_with_mysql import *


def calculate_average_lab_study(intdf,efcoldf):
    '''
    Computes the average emission factors (EF) for lab studies.
    This function computes the average EF for each fire type and study combination. It excludes specific fire types 
    such as 'cookstove', 'open cooking', etc. If there are multiple EF columns for a fire type-study pair, 
    an average column is added to intdf and the individual columns are dropped. The efcoldf dataframe is also 
    updated accordingly.

    Parameters:
    - intdf (pd.DataFrame): Input dataframe containing emission values.
    - efcoldf (pd.DataFrame): Information dataframe about EF columns.

    Returns:
    - pd.DataFrame: Updated intdf with averaged EF columns.
    - pd.DataFrame: Updated efcoldf with information on the new averaged columns.

    Note:
    Assumes a specific structure and content in the input dataframes.
    '''
    # getting the list of unique fire types from efcoldf(information table of all ef columns of pdb)
    efcoldf=sort_by_fire_type_year(efcoldf)
    fire_type_ll=efcoldf['fire_type'].unique()
    fire_type_ll = set(fire_type_ll) - {'cookstove', 'open cooking', 'dung burning', 'charcoal burning','coal burning'}
    for ft in fire_type_ll:
        study_ll=efcoldf['study'][efcoldf['fire_type']==ft].unique().tolist() # the list of unique study of a specific fire type.
        df=pd.DataFrame()
        for study in study_ll:
            efcols=efcoldf['efcol'][efcoldf['fire_type']==ft][efcoldf['study']==study][efcoldf['measurement_type']=='lab'].to_list()
            if len(efcols)>1:
                newcol='EF_'+ft.replace(' ','_')+'_'+study
                intdf[newcol]=intdf[efcols].T.mean()
                intdf=intdf.drop(columns=efcols)
                df['efcol']=[newcol]
                df['fire_type']=[ft]
                df['study']=[study]
                df['fuel_type']=','.join(efcoldf['fuel_type'][efcoldf['efcol'].isin(efcols)].unique())
                df['measurement_type']=['lab']                       
                efcoldf=efcoldf[~efcoldf['efcol'].isin(efcols)]
                efcoldf = pd.concat([efcoldf, df], ignore_index=True)
                efcoldf=efcoldf.reset_index(drop=True)
    return intdf, efcoldf


def get_lab_study_fuel_avg(intdf,efcoldf):
    '''
    Computes the average emission factors (EF) for lab studies based fuel types of specific fire type.

    This function calculates the average EF for each fire type and fuel type combination in lab studies. 
    If there are multiple EF columns for a study-fuel type pair for a specific fire type, the function computes an average column 
    and adds it to intdf, while the individual columns are removed. The efcoldf dataframe is also updated with 
    information on the new averaged columns.

    Parameters:
    - intdf (pd.DataFrame): Input dataframe containing emission values.
    - efcoldf (pd.DataFrame): Information dataframe about EF columns.

    Returns:
    - pd.DataFrame: Updated intdf with averaged EF columns for specific fuel types.
    - pd.DataFrame: Updated efcoldf with information on the new averaged columns.

    Note:
    Assumes a specific structure and content in the input dataframes.
    '''
    # Getting the list of unique fire types from efcoldf(information table of all ef columns of pdb)
    efcoldf=sort_by_fire_type_year(efcoldf)
    fire_type_ll=efcoldf['fire_type'].unique().tolist()
    for ft in fire_type_ll:
        fuel_ll=efcoldf['fuel_type'][efcoldf['fire_type']==ft][efcoldf['measurement_type']=='lab'].unique().tolist() # The list of unique fuel
        df=pd.DataFrame()
        for fuel in fuel_ll:
            efcols=efcoldf['efcol'][efcoldf['fire_type']==ft][efcoldf['fuel_type']==fuel][efcoldf['measurement_type']=='lab'].to_list()
            st_ll=list(efcoldf['study'][efcoldf['efcol'].isin(efcols)].unique())
            legend_ll=list(efcoldf['legend'][efcoldf['efcol'].isin(efcols)].unique())
            newcol='EF_'+fuel.replace(' ','_')+'_avg'
            
            intdf[newcol]=intdf[efcols].T.mean()
            intdf=intdf.drop(columns=efcols)
            
            df['efcol']=[newcol]
            df['fire_type']=[ft]
            df['study']=','.join(st_ll)
            df['measurement_type']=['lab']
            df['legend']=','.join(legend_ll)                  
            
            efcoldf=efcoldf[~efcoldf['efcol'].isin(efcols)]
            efcoldf = pd.concat([efcoldf, df], ignore_index=True)
            efcoldf=efcoldf.reset_index(drop=True)
    return intdf, efcoldf



def drop_cols_er_adj_calc(df, field, lab):
    '''
    Filters filed/lab list for emission ratio (ER) adjustment calculation.
    if a EF column with 'inorganic gas', 'methane', and 'NMOC_g' pollutant categories
    has all null-values than it is excluded from the input 'field'/'lab' list.

    Parameters:
    - df (pd.DataFrame): Input dataframe with pollutant data.
    - field (list): Columns representing field data.
    - lab (list): Columns representing lab data.
    '''
    
    gases=['inorganic gas', 'methane', 'NMOC_g']
    df_test=df[df['pollutant_category'].isin(gases)].reset_index(drop=True)
    
    new_field=[]
    for col in field:
        if len(df_test[df_test[col].notna()])!=0:
            new_field.append(col)
    
    new_lab=[]
    for col in lab:
        if len(df_test[df_test[col].notna()])!=0:
            if col!='EF_akagi11_africa_stockwell15': # Eliminating this lab data because it is already adjusted
                new_lab.append(col)
    return new_field, new_lab


def lab_data_adjust_to_field_conditions (df,efcoldf):
    '''
   Emission Ratio adjustment calculation.
    
    Parameters:
    - df: A DataFrame containing emission data.
    - efcoldf: A DataFrame containing information about emission factors (EF).

    Returns:
    - Modified 'df' with adjusted ER values.
    - Updated 'efcoldf' with new column names for ER-adjusted values.
    '''
    bk_db=connect_db('backend_db')
    # Iterate through different fire types excluding 'peat'.
    for fire_type in set(efcoldf['fire_type'].unique())-{'peat'}: 
        # Import a list compounds emitted during flaming combustion.
        fl=pd.read_sql(text('select * from bkdb_compound_flaming_combustion_type'), con=bk_db)
        # Drop BC from the list
        fl=fl.drop(index=8)
        
        # Extract the list of lab and field EF columns for a specific fire type.
        lab=efcoldf['efcol'][efcoldf['fire_type']==fire_type][efcoldf['measurement_type']=='lab'].tolist()
        field=efcoldf['efcol'][efcoldf['fire_type']==fire_type][efcoldf['measurement_type']=='field'].tolist()
        
        # Update lab and field Ef column list.
        field=drop_cols_er_adj_calc(df, field, lab)[0]
        lab=drop_cols_er_adj_calc(df, field, lab)[1]
        
        # Get inidces for CO and CO2.
        co_ind=get_ind(df)[0] 
        co2_ind=get_ind(df)[1] 
        
        # Eliminate lab col if it doesn't have CO
        ll=set(df[lab].iloc[co_ind].index)
        ll2=set(df[lab].iloc[co_ind].dropna().index)
        drop_col=ll-ll2
        lab=list(set(lab)-drop_col)

        
        # ER ADJ correction is applied to gases.
        gases=['inorganic gas', 'methane', 'NMOC_g']
        
        # Calculate field study averages
        fieldavg=df[field].T.mean()
        
        for l in lab:
            new_col=l+'_ER_ADJ'
            lab_ind=set(df[l][df['pollutant_category'].isin(gases)].dropna().index)
            
            fl_ind=set(df[df.index.isin(lab_ind)][df['id'].isin(fl['id'])].index)
            sml_ind=lab_ind-fl_ind
            
            co2_cor=df[l][df.index.isin(fl_ind)]/df[l].iloc[co2_ind]
            co2_cor[new_col]=co2_cor * fieldavg.iloc[co2_ind]
            
            co_cor=df[l][df.index.isin(sml_ind)]/df[l].iloc[co_ind]
            co_cor[new_col]=co_cor * fieldavg.iloc[co_ind]
            
            cor = pd.DataFrame(pd.concat([co2_cor[new_col], co_cor[new_col]]))
            cor.columns = [new_col]
            
            df=df.merge(cor, right_index=True, left_index=True, how='left')
            
            iind=list(df[df[l].notna()][df[new_col].isnull()].index)
            for xx in iind:
                df.loc[xx,new_col]=df[l].iloc[xx]
            # Changing col name in df and efcoldf    
            df=df.drop(columns=l)
            alt_ind=efcoldf[efcoldf['efcol']==l].index
            efcoldf.loc[alt_ind,'efcol']=new_col
    
    return df, efcoldf
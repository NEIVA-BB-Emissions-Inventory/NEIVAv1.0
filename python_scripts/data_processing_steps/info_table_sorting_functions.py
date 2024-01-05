#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 11:52:54 2023

@author: samiha
"""

import pandas as pd
from NEIVA.python_scripts.connect_with_mysql import connect_db
from sqlalchemy import text

def sort_by_fire_type_year(data):
    sorting_list=['savanna', 'boreal forest','tropical forest', 'temperate forest',
              'peatland','peat', 'chaparral', 'open cooking', 'cookstove', 'dung burning',
              'charcoal making', 'charcoal burning', 'coal burning','pasture maintenance', 
              'crop residue', 'garbage burning','grassland', 'multiple fire type']
    df=pd.DataFrame()
    
    for ft in sorting_list:
        if len(data.columns[data.columns.isin(['year_akagi_data'])])==1:
            dd=data[data['fire_type']==ft]
            dd=dd.sort_values(by=['year','year_akagi_data'], ascending=False)
        if len(data.columns[data.columns.isin(['year_akagi_data'])])!=1:
            dd=data[data['fire_type']==ft].sort_values(by='year', ascending=False)
        
        dd=dd.reset_index(drop=True)
        df=pd.concat([df,dd],ignore_index=True)
        df=df.reset_index(drop=True)
        
    return df

def sort_info_table_name (data):
    ldb=data[data['db']=='ldb'].reset_index(drop=True)
    rdb=data[data['db']=='rdb'].reset_index(drop=True)
    pdb=data[data['db']=='pdb'].reset_index(drop=True)
    
    ldb=sort_by_fire_type_year(ldb)
    rdb=sort_by_fire_type_year(rdb)
    pdb=sort_by_fire_type_year(pdb)
    
    data_final = pd.concat([ldb, rdb, pdb], ignore_index=True)
    
    return data_final
    


def assign_year_col_efcoldf(efcoldf):
    for i in efcoldf[~efcoldf['efcol'].str.contains('kagi')].index:
        efcoldf.loc[i,'year']=float(efcoldf['efcol'].iloc[i][-2:])
    
    for i in efcoldf[efcoldf['efcol'].str.contains('kagi')].index:
        efcoldf.loc[i,'year']=11
        efcoldf.loc[i,'year_akagi_data']=float(efcoldf['efcol'].iloc[i].replace(')','')[-2:])
        
    for i in range(len(efcoldf)):
        efcoldf.loc[i,'year']=float('20'+str(efcoldf['year'].iloc[i]))
    
    
    for i in efcoldf[efcoldf['year_akagi_data'].isin([3,0,8,7,6,9,4,1])].index:
        efcoldf.loc[i,'year_akagi_data']=float('200'+str(efcoldf['year_akagi_data'].iloc[i]))
    
    
    for i in efcoldf[efcoldf['year_akagi_data'].isin([11,13,15,10])].index:
        efcoldf.loc[i,'year_akagi_data']=float('20'+str(efcoldf['year_akagi_data'].iloc[i]))
    
    
    for i in efcoldf[efcoldf['year_akagi_data'].isin([97,93,91,95,96,98,94,99])].index:
        efcoldf.loc[i,'year_akagi_data']=float('19'+str(efcoldf['year_akagi_data'].iloc[i]))
    return efcoldf


def assign_legend_col(efcoldf):
    
    for i in efcoldf[efcoldf['efcol'].str.contains('kagi')].index:
        sst=efcoldf['study'].iloc[i].replace('Akagi_11(','').replace(')','')
        ll=sst[0].upper()+sst[1:-2]+' et al '+sst[-2:]
        efcoldf.loc[i,'legend']='Akagi_11('+ll+')'
    
    for i in efcoldf[~efcoldf['efcol'].str.contains('kagi')].index:
        sst=efcoldf['study'].iloc[i]
        ll=sst[0].upper()+sst[1:-2]+' et al '+sst[-2:]
        efcoldf.loc[i,'legend']=ll
        
    return efcoldf
    
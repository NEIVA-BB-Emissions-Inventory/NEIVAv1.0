#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 11:52:54 2023

@author: samiha
"""

import pandas as pd
from connect_with_mysql import*

con=connect_db('NEIVA_db')


data=pd.read_sql('select * from info_table_name', con=con)
efcoldf=pd.read_sql('select * from info_table', con=con)

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
        df=df.append(dd)
        df=df.reset_index(drop=True)
        
    return df

def sort_info_table_name (data):
    ldb=data[data['db']=='ldb'].reset_index(drop=True)
    rdb=data[data['db']=='rdb'].reset_index(drop=True)
    pdb=data[data['db']=='pdb'].reset_index(drop=True)
    
    ldb=sort_by_fire_type_year(ldb)
    rdb=sort_by_fire_type_year(rdb)
    pdb=sort_by_fire_type_year(pdb)
    
    data_final=ldb.append(rdb).append(pdb)
    
    return data_final
    

def assign_multiple_fire_type (data, efcoldf):
    tt=pd.DataFrame()
    
    pdb=data[data['db']=='pdb']
    pdb_tbl=pdb['tbl_name'][pdb['fire_type']=='multiple fire type'].to_list()
    
    rdb_tbl=[]
    for tbl in pdb_tbl:
        rdb_tbl.append(tbl.replace('pdb','rdb'))
    
    all_tbl=pdb_tbl+rdb_tbl
    
    tt['tbl_name']=all_tbl
    
    for i in tt[tt['tbl_name'].str.contains('pdb')].index:
        tbl=tt['tbl_name'].iloc[i]
        dat=pd.read_sql('select * from '+tbl, con=con)
        efcol=list(dat.columns[dat.columns.str.contains('EF')])
        
        ft_ll=list(efcoldf['fire_type'][efcoldf['efcol'].isin(efcol)].unique())
        tt.loc[i,'fire_type']=', '.join(ft_ll)
        r_ind=tt[tt['tbl_name']==tbl.replace('pdb','rdb')].index[0]
        tt.loc[r_ind,'fire_type']=', '.join(ft_ll)

    data=data.merge(tt[['tbl_name','fire_type']], on='tbl_name', how='left')   
    
    ii=list(data[data['fire_type_y'].notna()].index)
    
    for i in ii:
        data.loc[i,'fire_type_x']=data['fire_type_y'].iloc[i]
    
    data['fire_type']=data['fire_type_x']
    data=data.drop(columns=['fire_type_x','fire_type_y'])    
    return data


def assign_pollutant_catagory (data):
    
    for i in range(len(data)):
        tbl_name=data['tbl_name'].iloc[i]
        print(tbl_name)
        dd=pd.read_sql('select * from '+tbl_name, con=con)
        dd=dd[dd['pollutant_category'].notna()].reset_index(drop=True)
        
        ind=dd[dd['pollutant_category'].str.contains('NMOC_p')].index
        for k in ind:
            dd.loc[k,'pollutant_category']='NMOC_p'
        
        ind=dd[dd['pollutant_category'].str.contains('PM')].index
        for k in ind:
            dd.loc[k,'pollutant_category']='PM'

        pl=list(dd['pollutant_category'][dd['pollutant_category'].notna()].unique())
        
        reference_list = ['inorganic gas','methane', 'NMOC_g', 'PM','NMOC_p']
        sorted_list = sorted(pl, key=lambda x: reference_list.index(x))

        print(sorted_list)
        
        data.loc[i,'pollutant category']=', '.join(sorted_list)
        
    return data
    

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













    
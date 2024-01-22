#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 10:25:48 2024

@author: samiha
"""

import pandas as pd
import numpy as np
import pubchempy as pcp
import matplotlib.pyplot as plt
from connect_with_mysql import*


bk_db=connect_db('backend_db')
primary_db=connect_db('primary_db')
legacy_db=connect_db('legacy_db')

dd=pd.read_sql('select * from bkdb_info_table_name', con=bk_db)

ind=list(dd['tbl_name'][dd['db']=='ldb'].index)

dd.loc[10,'tbl_name']='ldb_peat'

for i in ind:
    print(dd['tbl_name'].iloc[i])
    df=pd.read_sql('select * from '+dd['tbl_name'].iloc[i], con=legacy_db)
    ll=list(df['pollutant_category'][df['pollutant_category'].notna()].unique())
    print(','.join(ll))
    dd.loc[i,'pollutant_category']=','.join(ll)


ind=list(dd['tbl_name'][dd['db']=='pdb'].index)

for i in ind:
    print(dd['tbl_name'].iloc[i])
    df=pd.read_sql('select * from '+dd['tbl_name'].iloc[i], con=primary_db)
    ll=list(df['pollutant_category'][df['pollutant_category'].notna()].unique())
    print(','.join(ll))
    dd.loc[i,'pollutant_category']=','.join(ll)
    
    
    print(dd['tbl_name'].iloc[i].replace('pdb','rdb'))
    ii=dd[dd['tbl_name']==dd['tbl_name'].iloc[i].replace('pdb','rdb')].index
    if len(ii)==1:
        dd.loc[ii[0],'pollutant_category']=','.join(ll)
   

dd.to_sql(name='bkdb_info_table_name',con=bk_db, if_exists='replace', index=False)

#___
pp=pd.read_sql('select * from property_surrogate_info', con=bk_db)

output_db=connect_db('neiva_output_db')
dd=pd.read_sql('select * from Property_Surrogate', con=output_db)


ii=list(dd[['compound','MOZT1']][dd['S07']=='TERP'][dd['S18B']=='SESQ'].index)

for i in ii:
    dd.loc[i,'MOZT1']='BCARY'



ii=list(dd[['compound','MOZT1']][dd['S07']=='TERP'][dd['MOZT1'].isnull()].index)

for i in ii:
    dd.loc[i,'MOZT1']='MTERP'



dd.to_sql(name='Property_Surrogate',con=output_db, if_exists='replace', index=False)





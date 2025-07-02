f#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 13:33:41 2024

@author: samiha
"""


import pubchempy as pcp
import pandas as pd
import numpy as np
from connect_with_mysql import connect_db

bk_db=connect_db('backend_db')
output_db=connect_db('neiva_output_db')

def get_ind(df, compound):
    if compound == 'PM2.5*':
        ind = df[df['compound'].str.contains('PM')][df['id']!='PM10'][df['id']!='PM2.5_ipcc'].index.tolist()
    elif compound in ['PM10', 'PM2.5','PM1','PM2.5(PM1-PM5)','OA', 'EC', 'BC', 'OC', 'NOx_as_NO']:
        ind = df[df['id'] == compound].index.tolist()
    else:
        iid = pcp.get_compounds(compound, 'name')[0].inchi
        ind = df[df['id'] == iid].index.tolist()
    return ind


df=pd.read_sql('select * from Processed_EF', con=output_db)
efcoldf=pd.read_sql('select * from info_efcol_processed_data', con=bk_db)
#allcol= ['legend','fuel_type','measurement_type','MCE','EF']
com_name='sulphur dioxide'
ind = get_ind (df, com_name)

# ind=df[df['id']=='NOx_as_NO'].index
# ind=[5]

efcol=list(efcoldf['efcol'])
efcoldf['EF']=df[efcol][df.index.isin(ind)].mean().values
ll=efcoldf[efcoldf['EF'].notna()]
ll=ll.reset_index(drop=True)

ll_final=pd.DataFrame()
ft_ll=['savanna', 'boreal forest', 'tropical forest', 'temperate forest','chaparral', 'peat', 'crop residue', 'garbage burning']
for i in range(len(ft_ll)):
    ldf=pd.DataFrame()
    ldf=ll[ll['fire_type']==ft_ll[i]]
    ldf=ldf.sort_values(by='EF', ascending=False)
    ll_final=ll_final.append(ldf)
    
#ll_final=ll_final.applymap(lambda x: rounding(x))
dd=ll_final[['legend', 'fuel_type', 'measurement_type','fire_type', 'EF']]
dd.to_excel('/Users/samiha/Desktop/send_data_kelley/SO2.xlsx',index=False)



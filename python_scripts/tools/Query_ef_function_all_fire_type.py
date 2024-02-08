#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 13:33:41 2024

@author: samiha
"""

bk_db=connect_db('backend_db')
output_db=connect_db('neiva_output_db')

df=pd.read_sql('select * from Processed_EF', con=output_db)
efcoldf=pd.read_sql('select * from info_efcol_processed_data', con=bk_db)
#allcol= ['legend','fuel_type','measurement_type','MCE','EF']
com_name='BC'
ind = get_ind (df, com_name)
efcol=list(efcoldf['efcol'])
efcoldf['EF']=df[efcol][df.index.isin(ind)].mean().values
ll=efcoldf[efcoldf['EF'].notna()]
ll=ll.reset_index(drop=True)

ll_final=pd.DataFrame()
ft_ll=['tropical forest', 'temperate forest', 'boreal forest', 'savanna', 'crop residue', 'peat']
for i in range(len(ft_ll)):
    ldf=pd.DataFrame()
    ldf=ll[ll['fire_type']==ft_ll[i]]
    ldf=ldf.sort_values(by='EF', ascending=False)
    ll_final=ll_final.append(ldf)
    
ll_final=ll_final.applymap(lambda x: rounding(x))
dd=ll_final[['legend', 'fuel_type', 'measurement_type','fire_type', 'EF']]
dd.to_excel('/Users/samiha/Desktop/send_data_kelley/BC.xlsx',index=False)


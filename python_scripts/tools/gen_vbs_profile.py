#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 19:30:13 2023

@author: samiha
"""
import pubchempy as pcp
import pandas as pd
import numpy as np
from NEIVA.python_scripts.connect_with_mysql import*
from NEIVA.python_scripts.tools.assign_mozart_species import mozart_species
from sqlalchemy import text
import matplotlib.pyplot as plt


def calc_VBS (dd, ft):
    output_db=connect_db('neiva_output_db')
    bk_db=connect_db('backend_db')
    
    nmog=dd[dd['pollutant_category']=='NMOC_g']
    
    #__
    unk1=nmog[nmog['compound'].str.contains('unk')]
    unk2=nmog[nmog['compound'].str.contains('Unk')]
    
    # Compounds with 'unknown' (unknown chemical structure)
    nmog=nmog[~nmog['compound'].str.contains('unk')]
    nmog=nmog[~nmog['compound'].str.contains('Unk')]
    nmog=nmog.reset_index(drop=True)
    
    # Loading dataset with chemmical mechanism and property assignments
    tt=pd.read_sql(text('select * from chem_property_inchi'), con=bk_db)
    tt1=pd.read_sql(text('select * from chem_property_h15isomers'), con=bk_db)
    tt2=pd.read_sql(text('select * from chem_property_lumpCom'), con=bk_db)
    
    #__
    tt2['cstar']=tt2[['cstar','cstar1','cstar2']].T.mean()
    #__
    
    # Final chem_property datasets
    tt_f=tt.append(tt1).append(tt2)
    tt_f=tt_f.reset_index(drop=True)
    tt_f=tt_f[['id','S07','S18B','S07T','MOZT1','S22', 'cstar']]
    
    nmog=nmog.merge(tt_f, on='id', how='left')
    nmog=mozart_species(nmog)
    
    ft=ft.replace(' ','_')
    efcol='AVG_'+ft
    nmog['ef']=nmog[efcol]
    
    nmog=nmog[nmog['ef'].notna()].reset_index(drop=True)
    
    # Load Lumped Compounds with Speciation dataset.
    lc_spec_ref= pd.read_sql(text('select * from chem_property_lumpCom_spec'), con=bk_db)
    # Get the Lumped compound id (ids without InChI)
    lc_spec_ref_id=list(lc_spec_ref['id'][~lc_spec_ref['id'].str.contains('InChI')])
    # Get the dataframe where lumped compound ids of 'chem_property_lumpCom_spec' is in 'nmogdf'
    lc_spec_df=nmog[nmog['id'].isin(lc_spec_ref_id)] 
    lc_spec_df=lc_spec_df.reset_index(drop=True)
    
    # Exclude the lumped compounds that has speciation from the 'nmogdf' dataset.
    nmog_final=nmog[~nmog['id'].isin(lc_spec_df['id'])]
    nmog_final=nmog_final.reset_index(drop=True)
    
    nmog_final=nmog_final[['ef', 'cstar']]
    #__
    df_final=pd.DataFrame() # The dataframe where the assignments of lumped compound from speceation will be assigned.
    # iterate over 'lc_spec_df' dataframe
    for i in range(len(lc_spec_df)):
        df=pd.DataFrame()
        # Search formulas of 'lc_spec_df' id 'lc_spec_ref'
        ll=lc_spec_ref[lc_spec_ref['formula']==lc_spec_df['formula'][i]]
        # Exclud the lumped compound id
        ll=ll[ll['id']!=lc_spec_df['id'][i]]
        ll=ll.reset_index(drop=True)
        # Iterate over the speciation dataset 
        for k in range(len(ll)):
            # Equally dividing the EF over the number of specie
            df.loc[k,'ef']=lc_spec_df['ef'].iloc[i]/len(ll)
            df.loc[k,'cstar']=ll['cstar'][k]
        df_final=df_final.append(df)
        df_final=df_final.reset_index(drop=True)
    # Adding the two datasets        
    nmog_final=nmog_final.append(df_final)        
    nmog_final=nmog_final.reset_index(drop=True)
    
    unk_totef=unk1[efcol].sum()+unk2[efcol].sum()
    nmog_final['ef']=nmog_final['ef']+(unk_totef/len(nmog_final))
    
    nmog_final['bin']=pd.cut(nmog_final['cstar'], bins=range(0,14,1))
    aa=nmog_final.groupby('bin')['ef'].sum()/nmog_final['ef'].sum()

    return aa





#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 21:00:30 2023

@author: samiha
"""

import pubchempy as pcp
import pandas as pd
import numpy as np
from NEIVA.python_scripts.connect_with_mysql import*
from NEIVA.python_scripts.tools.assign_mozart_species import mozart_species
from sqlalchemy import text
import matplotlib.pyplot as plt


def calc_OHR (dd,chem, tot_voc, ft):
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
    
    # Final chem_property datasets
    tt_f=tt.append(tt1).append(tt2)
    tt_f=tt_f.reset_index(drop=True)
    tt_f=tt_f[['id','S07','S18B','S07T','MOZT1','S22', 'kOH']]
    
    nmog=nmog.merge(tt_f, on='id', how='left')
    nmog=mozart_species(nmog)

    ft=ft.replace(' ','_')
    efcol='AVG_'+ft
    nmog['ef']=nmog[efcol]
    nmog['mole']=nmog['ef']/nmog['mm']
    
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
    
    nmog_final=nmog_final[[chem,'mole']]
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
            df.loc[k,chem]=ll[chem][k]
            # Equally dividing the EF over the number of specie
            df.loc[k,'mole']=lc_spec_df['mole'].iloc[i]/len(ll)
            #df.loc[k,'kOH']=ll['kOH'][k]
        df_final=df_final.append(df)
        df_final=df_final.reset_index(drop=True)
    # Adding the two datasets        
    nmog_final=nmog_final.append(df_final)        
    nmog_final=nmog_final.reset_index(drop=True)
    
    unk_tot_mole=(unk1[efcol]/unk1['mm']).sum()+(unk2[efcol]/unk2['mm']).sum()
    
    nmog_final['mole']=nmog_final['mole']+(unk_tot_mole/len(nmog_final))
    
    nmog_final['mole_frac']=nmog_final['mole']/nmog_final['mole'].sum()
    nmog_final['conc']=nmog_final['mole_frac']*tot_voc # tot_voc in ppb
    nmog_final['ohr']=nmog_final['conc']*2.5e10*nmog_final['kOH']      

    uu=list(nmog_final[chem][nmog_final[chem].notna()].unique())
    # Creating the model species profile dataframe
    ohrdf=pd.DataFrame()
    for i in range(len(uu)):
        ohrdf.loc[i,chem]=uu[i]
        ohrdf.loc[i,'OHR']=nmog_final['ohr'][nmog_final[chem]==uu[i]].sum()
    
    tot_ohr=ohrdf['OHR'].sum()
    ohrdf['OHR_frac']=ohrdf['OHR']/tot_ohr
    
    ohrdf=ohrdf.sort_values(by='OHR', ascending=False)
    ohrdf=ohrdf.reset_index(drop=True)

    return ohrdf

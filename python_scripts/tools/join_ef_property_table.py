#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 10:18:51 2024

@author: samiha
"""

import pandas as pd
import numpy as np
import pubchempy as pcp
import matplotlib.pyplot as plt
from NEIVA.python_scripts.connect_with_mysql import*
from NEIVA.python_scripts.tools.assign_mozart_species import mozart_species
from NEIVA.python_scripts.tools.assign_geos_chem_species import geos_chem_species
from NEIVA.python_scripts.data_integration_process.sort_molec_formula import *
from sqlalchemy import text

def join_ef_property (dd):
    output_db=connect_db('neiva_output_db')
    bk_db=connect_db('backend_db')
    
    nmog=dd[dd['pollutant_category']=='NMOC_g']
    nmog=nmog[~nmog['compound'].str.contains('unk')]
    nmog=nmog[~nmog['compound'].str.contains('Unk')]
    nmog=nmog.reset_index(drop=True)
    # Loading dataset with chemical mechanism and property assignments
    tt=pd.read_sql(text('select * from chem_property_inchi'), con=bk_db)
    tt1=pd.read_sql(text('select * from chem_property_h15isomers'), con=bk_db)
    tt2=pd.read_sql(text('select * from chem_property_lumpCom'), con=bk_db)
    # Final chem_property datasets
    tt_f = pd.concat([tt, tt1, tt2], ignore_index=True)
    tt_f=tt_f.reset_index(drop=True)
    tt_f=tt_f[['id','S07','S18B','S07T','MOZT1','S22','GEOS_chem',\
               'vp_nannoolal','vp','cstar','kOH','ko3_exp',\
                'kno3_exp','hc']]
    
    nmog=nmog.merge(tt_f, on='id', how='left')
    nmog=mozart_species(nmog)
    nmog=geos_chem_species(nmog)
   
    return nmog

def assign_GFED_lumed_ef_conv_factor (nmog):
    
    for i in range(len(nmog)):
         nmog.loc[i,'nC']=get_nMolecule('C',nmog['formula'].iloc[i])
    nmog['conversion_factor']= (nmog['nC']*12)/nmog['mm']
    
    return nmog

def lump_com_with_speciation (nmog, chem):
    bk_db=connect_db('backend_db')
    # Assign GFED cinversion factor
    nmog = assign_GFED_lumed_ef_conv_factor (nmog)
    # Load Lumped Compounds with Speciation dataset.
    lc_spec_ref= pd.read_sql(text('select * from chem_property_lumpCom_spec'), con=bk_db)
    lc_spec_ref = assign_GFED_lumed_ef_conv_factor (lc_spec_ref)
    # Get the Lumped compound id which is id's without an InChI
    lc_spec_ref_id=list(lc_spec_ref['id'][~lc_spec_ref['id'].str.contains('InChI')])
    # Get the dataframe where lumped compound ids of 'chem_property_lumpCom_spec' is in 'nmogdf'
    lc_spec_df=nmog[nmog['id'].isin(lc_spec_ref_id)] 
    lc_spec_df=lc_spec_df.reset_index(drop=True)
    
    # Exclude the lumped compounds that has speciation from the 'nmogdf' dataset.
    nmog_final=nmog[~nmog['id'].isin(lc_spec_df['id'])]
    nmog_final=nmog_final.reset_index(drop=True)

    nmog_final=nmog_final[['mm',chem,'ef','conversion_factor','kOH','cstar','vp_nannoolal','vp','hc','ko3_exp','kno3_exp']]
    
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
            df.loc[k,'ef']=lc_spec_df['ef'].iloc[i]/len(ll)
            df.loc[k,'mm']=ll['mm'][k]
            df.loc[k,'conversion_factor']=ll['conversion_factor'][k]
            df.loc[k,'kOH']=ll['kOH'][k]
            df.loc[k,'cstar']=ll['cstar'][k]
            df.loc[k,'vp_nannoolal']=ll['vp_nannoolal'][k]
            df.loc[k,'vp']=ll['vp'][k]
            df.loc[k,'hc']=ll['hc'][k]
            df.loc[k,'ko3_exp']=ll['ko3_exp'][k]
            df.loc[k,'kno3_exp']=ll['kno3_exp'][k]
        df_final = pd.concat([df_final, df], ignore_index=True)
        df_final=df_final.reset_index(drop=True)
    # Adding the two datasets        
    nmog_final = pd.concat([nmog_final, df_final], ignore_index=True)        
    nmog_final=nmog_final.reset_index(drop=True)
    
    return nmog_final


def distribute_unk_ef (dd, efcol, nmog_final):
    
    nmog=dd[dd['pollutant_category']=='NMOC_g']
    unk1=nmog[nmog['compound'].str.contains('unk')]
    unk2=nmog[nmog['compound'].str.contains('Unk')]
    
    tot_unk_ef = unk1[efcol].sum() + unk2[efcol].sum()
    sum_factor = tot_unk_ef/len(nmog_final)
    
    nmog_final['ef']=nmog_final['ef']+ sum_factor
    
    return nmog_final

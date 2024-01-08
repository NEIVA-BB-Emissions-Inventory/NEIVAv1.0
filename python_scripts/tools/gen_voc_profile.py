#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 11:53:38 2023

@author: samiha
"""

import pubchempy as pcp
import pandas as pd
import numpy as np
from NEIVA.python_scripts.connect_with_mysql import*
from NEIVA.python_scripts.tools.assign_mozart_species import mozart_species
from sqlalchemy import text
import matplotlib.pyplot as plt

            
def voc_profile(dd, chem, ft):
    '''
    Calculates the VOC for a specific chemcial mechanism and fire type.
    '''
    output_db=connect_db('neiva_output_db')
    bk_db=connect_db('backend_db')
    
    # Loading Recommended_EF data, extracting NMOG_C and excluding compounds with 'unknown' (unknown chemical structure)
    nmog=dd[dd['pollutant_category']=='NMOC_g']
    
    unk1=nmog[nmog['compound'].str.contains('unk')]
    unk2=nmog[nmog['compound'].str.contains('Unk')]
    
    nmog=nmog[~nmog['compound'].str.contains('unk')]
    nmog=nmog[~nmog['compound'].str.contains('Unk')]
    nmog=nmog.reset_index(drop=True)
    
    # Loading dataset with chemical mechanism and property assignments
    tt=pd.read_sql(text('select * from chem_property_inchi'), con=bk_db)
    tt1=pd.read_sql(text('select * from chem_property_h15isomers'), con=bk_db)
    tt2=pd.read_sql(text('select * from chem_property_lumpCom'), con=bk_db)
    
    # Checking if all ids of NMOC_gs is matched with the chem_property datasets or not
    # id_all=list(tt['id'])+list(tt1['id'])+list(tt2['id'])+list(tt3['id'])
    # assert len(nmog[['mm','formula','compound','id']][~nmog['id'].isin(id_all)])==0
    
    # Final chem_property datasets
    tt_f = pd.concat([tt, tt1, tt2], ignore_index=True)
    tt_f=tt_f.reset_index(drop=True)
    tt_f=tt_f[['id','S07','S18B','S07T','MOZT1','S22']]
    
    nmog=nmog.merge(tt_f, on='id', how='left')
    
    nmog=mozart_species(nmog)
    
    # Set EF column based on the input parameter 'fire type'
    efcol='AVG_'+ft.replace(' ','_')
    nmog['ef']=nmog[efcol]
    
    nmog=nmog[nmog['ef'].notna()].reset_index(drop=True)

    # Load Lumped Compounds with Speciation dataset.
    lc_spec_ref= pd.read_sql(text('select * from chem_property_lumpCom_spec'), con=bk_db)
    # Get the Lumped compound id which is id's without an InChI
    lc_spec_ref_id=list(lc_spec_ref['id'][~lc_spec_ref['id'].str.contains('InChI')])
    # Get the dataframe where lumped compound ids of 'chem_property_lumpCom_spec' is in 'nmogdf'
    lc_spec_df=nmog[nmog['id'].isin(lc_spec_ref_id)] 
    lc_spec_df=lc_spec_df.reset_index(drop=True)
    
    # Exclude the lumped compounds that has speciation from the 'nmogdf' dataset.
    nmog_final=nmog[~nmog['id'].isin(lc_spec_df['id'])]
    nmog_final=nmog_final.reset_index(drop=True)
    
    # Process the 'TERP' if chem is MOZART
    if chem=='MOZT1':
        terp_unassigned=nmog_final['ef'][nmog_final['S07']=='TERP'][nmog_final['MOZT1'].isnull()].sum()
    
    nmog_final=nmog_final[[chem,'ef','mm']]            
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
            df.loc[k,'ef']=lc_spec_df['ef'].iloc[i]/len(ll)
            df.loc[k,'mm']=ll['mm'][k]
        df_final = pd.concat([df_final, df], ignore_index=True)
        df_final=df_final.reset_index(drop=True)
    # Adding the two datasets        
    nmog_final = pd.concat([nmog_final, df_final], ignore_index=True)        
    nmog_final = nmog_final.reset_index(drop=True)
    #___
    # The list of unique model species
    uu=list(nmog_final[chem][nmog_final[chem].notna()].unique())
    
    # Creating the model species profile dataframe
    prdf=pd.DataFrame()
    for i in range(len(uu)):
        prdf.loc[i,chem]=uu[i]
        totef=nmog_final['ef'][nmog_final[chem]==uu[i]].sum()
        mean_mm=nmog_final['mm'][nmog_final[chem]==uu[i]].mean()
        prdf.loc[i,'ef']=totef
        prdf.loc[i,'mm']=mean_mm
    
    # The TERP total EF is distributed equally to the following model species.
    if chem=='MOZT1':
        m_terp_sp=['MYRC','APIN','BPIN','LIMON']
        for i in range(len(m_terp_sp)):
            ind=prdf[prdf['MOZT1']==m_terp_sp[i]].index[0]
            prdf.loc[ind,'ef']=prdf['ef'].iloc[ind]+terp_unassigned/len(m_terp_sp)
            
    prdf['mole']=prdf['ef']/prdf['mm']
    
    all_unk = (unk1[efcol]/unk1['mm']).sum()+ (unk2[efcol]/unk2['mm']).sum()
    
    prdf['mole']=prdf['mole']+(all_unk/len(prdf))
        
    totmole=prdf['mole'].sum()
    prdf['mole_fraction']=prdf['mole']/totmole
    
    prdf=prdf.sort_values(by='mole_fraction', ascending=False)
    prdf=prdf.reset_index(drop=True)
    
    prdf['mole']=round(prdf['mole'],4)    
    prdf['mole_fraction']=round(prdf['mole_fraction'],4)

    return prdf


def plot_voc_profile(pp):
    
    x=np.arange(len(pp))
    plt.barplot(x, pp['mole_fraction'])
    plt.xticks()
    
    
    return


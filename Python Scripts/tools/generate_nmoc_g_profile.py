#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 11:53:38 2023

@author: samiha
"""

import pubchempy as pcp
import pandas as pd
import numpy as np
from connect_with_mysql import*

import matplotlib.pyplot as plt
output_db=connect_db('neiva_output_db')
bk_db=connect_db('backend_db')

# Loading Recommended_EF data, extracting NMOG_C and excluding compounds with 'unknown' (unknown chemical structure)
dd=pd.read_sql('select * from Recommended_EF', con=output_db)
nmog=dd[dd['pollutant_category']=='NMOC_g']
nmog=nmog[~nmog['compound'].str.contains('unk')]
nmog=nmog[~nmog['compound'].str.contains('Unk')]
nmog=nmog.reset_index(drop=True)

# Loading dataset with chemmical mechanism and property assignments
tt=pd.read_sql('select * from chem_property_inchi', con=bk_db)
tt1=pd.read_sql('select * from chem_property_h15isomers', con=bk_db)
tt2=pd.read_sql('select * from chem_property_lumpCom', con=bk_db)

# Checking if all ids of NMOC_gs is matched with the chem_property datasets or not
# id_all=list(tt['id'])+list(tt1['id'])+list(tt2['id'])+list(tt3['id'])
# assert len(nmog[['mm','formula','compound','id']][~nmog['id'].isin(id_all)])==0

# Final chem_property datasets
tt_f=tt.append(tt1).append(tt2)
tt_f=tt_f.reset_index(drop=True)
tt_f=tt_f[['id','S07','S18B','S07T','MOZT1','S22']]

nmog=nmog.merge(tt_f, on='id', how='left')

# Unassigned MOZT1 and mapping S18B species to MOZT1 species
nmog['MOZT1'].unique()
nmog['S18B'][nmog['MOZT1'].isnull()].unique()

# S18 model species of unassigned compounds
ll=['HCHO','KET2','OTH3','AMINS','ACYLS','OTH4','NROG',
    'OLEP','ROOH','OLE1', 'OLE2','FURNS','RCOOH','RCHO',
    'R1NO3','LVKS','OLEC','ARO1','OLED','BUT13','OLEA2',
    'OLEA1','AFG1','STYRS','PHEN','OLE4','ALK5','MEK','CRES',
    'ARO2','OTH2','NAPS','DLIMO','SVPHE','MVK','XYNL']

# Corresponding MOZT1 model species
ll2=['CH2O','MEK','BIGALK','C3H6','BIGENE','BIGENE','NROG',
     'BIGENE','XYLENES', 'BIGENE', 'BIGENE','XYLENES','TOLUENE',
     'BIGALK','ALKNIT','MEK','BIGENE','XYLENES','BIGENE','BIGENE','MEK',
     'TOLUENE','MEK','BIGENE','PHENOL','BIGENE','BIGALK','MEK','CRESOL',
     'XYLENES','NROG','XYLENES','LIMON','PHENOL','MVK','CRESOL']

adf=pd.DataFrame()
adf['MOZRT_sp']=ll2
adf['S18B_sp']=ll

# Assign MOZT1 model species
for aa in range(len(adf)):
    s=adf['S18B_sp'].iloc[aa]
    m=adf['MOZRT_sp'].iloc[aa]
    ind=nmog[nmog['S18B']==s][nmog['MOZT1'].isnull()].index
    if len(ind)!=0:
        for ii in ind:
            nmog.loc[ii,'MOZT1']=m
            
def voc_profile(nmog, chem, ft):
    '''
    Calculates the VOC for a specific chemcial mechanism and fire type.
    '''
    # Set EF column based on the input parameter 'fire type'
    efcol='AVG_'+ft.replace(' ','_')
    nmog['ef']=nmog[efcol]
    
    nmog=nmog[nmog['ef'].notna()].reset_index(drop=True)

    # Load Lumped Compounds with Speciation dataset.
    lc_spec_ref= pd.read_sql('select * from chem_property_lumpCom_spec', con=bk_db)
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
        df_final=df_final.append(df)
        df_final=df_final.reset_index(drop=True)
    # Adding the two datasets        
    nmog_final=nmog_final.append(df_final)        
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
            
    totmole=prdf['mole'].sum()
    prdf['mole_fraction']=prdf['mole']/totmole
    
    prdf=prdf.sort_values(by='mole_fraction', ascending=False)
    prdf=prdf.reset_index(drop=True)
    
    prdf['mole']=round(prdf['mole'],4)    
    prdf['mole_fraction']=round(prdf['mole_fraction'],4)

    return prdf

chem='S07'
ft='temperate forest'
pp=voc_profile(nmog, chem, ft)    

#Save pathway
pathway='/Users/samiha/Desktop/NEIVA_final_v2/plot_data/neiva_temperate_forest_S07.xlsx'

pp.to_excel('/Users/samiha/Desktop/NEIVA_final_v2/plot_data/neiva_temperate_forest_S07.xlsx')

#pp.to_excel(pathway+'neiva_'+ft.replace(' ','_')+'_'+chem+'.xlsx',index=False)



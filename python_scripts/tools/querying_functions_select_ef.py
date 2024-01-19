#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 12:34:51 2024

@author: samiha
"""

import pandas as pd
import numpy as np
import pubchempy as pcp
import matplotlib.pyplot as plt
from NEIVA.python_scripts.connect_with_mysql import*
from sqlalchemy import text


# This function returns the PM2.5, OC, BC data for a specified fire type and table name. The table name inlcudes three different tables
# integrated ef, processed ef, recommended ef.
def select_pm_data (ft, table_name):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')
    if table_name=='integrated ef':
        efcoldf=pd.read_sql(text('select * from bkdb_info_efcol'), con=bk_db)
        df=pd.read_sql(text('select * from Integrated_EF'), con=output_db)

        efcol=list(efcoldf['efcol'])
        efcoldf['BC']=df[efcol][df['id']=='BC'].values[0]
        efcoldf['OC']=df[efcol][df['id']=='OC'].values[0]
        efcoldf['PM2.5']=df[efcol][df['id']=='PM2.5'].values[0]
        efcoldf['OA']=df[efcol][df['id']=='OC'].values[0]

        return (efcoldf[['legend','MCE','PM2.5','OC','BC','OA']][efcoldf['fire_type']==ft].reset_index(drop=True))

    if table_name=='processed ef':
        efcoldf=pd.read_sql(text('select * from info_efcol_processed_data'), con=bk_db)
        df=pd.read_sql(text('select * from Processed_EF'), con=output_db)

        efcol=list(efcoldf['efcol'])
        efcoldf['BC']=df[efcol][df['id']=='BC'].values[0]
        efcoldf['OC']=df[efcol][df['id']=='OC'].values[0]
        efcoldf['PM2.5']=df[efcol][df['id']=='PM2.5'].values[0]
        efcoldf['OA']=df[efcol][df['id']=='OC'].values[0]

        return (efcoldf[['legend','PM2.5','OC','BC','OA']][efcoldf['fire_type']==ft].reset_index(drop=True))

    if table_name=='recommended ef':
        df=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
        efcol='AVG_'+ft.replace(' ','_')
        iid=['PM<2.5','OC','BC','OA']
        return df[['compound',efcol]][df['id'].isin(iid)].reset_index(drop=True)


# This function returns the EF of specified pollutant category and fire type.
def select_ef_pollutant_category(ft, pc):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')
    
    efcol='AVG_'+ft.replace(' ','_')
    rdf=pd.read_sql(text('select * from Recommended_EF'), con=output_db)

    if pc=='PM optical property':
      return (rdf[['compound',efcol]][rdf['pollutant_category']==pc][rdf[efcol].notna()].reset_index(drop=True))
    else:
      return (rdf[['mm','formula','compound',efcol]][rdf['pollutant_category']==pc][rdf[efcol].notna()].reset_index(drop=True))

# This function returns the EF of specified compound name and table name (integrated ef, processed ef and recommended ef)
def select_compound(ft, com_name,table_name):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')

    if table_name=='integrated ef':
        df=pd.read_sql(text('select * from Integrated_EF'), con=output_db)
        efcoldf=pd.read_sql(text('select * from bkdb_info_efcol'), con=bk_db)
        allcol= ['legend','fuel_type','measurement_type','MCE',com_name]

    if table_name=='processed ef':
        df=pd.read_sql(text('select * from Processed_EF'), con=output_db)
        efcoldf=pd.read_sql(text('select * from info_efcol_processed_data'), con=bk_db)
        allcol= ['legend','measurement_type',com_name]

    aa=pcp.get_compounds(com_name, 'name')[0].inchi
    ind=df[df['id']==aa].index[0]

    efcol=list(efcoldf['efcol'])

    efcoldf[com_name]=df[efcol].iloc[ind].values

    ll=efcoldf[allcol][efcoldf['fire_type']==ft]
    ll=ll.sort_values(by='measurement_type')
    ll=ll[ll[com_name].notna()]
    ll=ll.reset_index(drop=True)
    return ll

def select_compound_rdf (ft, com_name):
    output_db=connect_db('neiva_output_db')
    df=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
    
    aa=pcp.get_compounds(com_name, 'name')[0].inchi
    ind=df[df['id']==aa].index[0]
    
    col='AVG_'+ft.replace(' ','_')
    
    return df[['mm','formula','compound',col]].iloc[ind]



# This function returns the speciation compounds of specified chemical mechanism and model surrogate
def speciation_profile(ft,chem,spc):
  bk_db=connect_db('backend_db')
  output_db=connect_db('neiva_output_db')
    
  pp=pd.read_sql(text('select * from Property_Surrogate'), con=output_db)
  df=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
  df=df.merge(pp[pp.columns[3:]],on='id', how='left')
  efcol='AVG_'+ft.replace(' ','_')
  return df[['formula','smile',efcol,chem,'cstar','kOH']][df[chem]==spc][df[efcol].notna()]

# Plots ef data of a specified fire type and table name
def plot_ef(compound,ft, table_name):
  bk_db=connect_db('backend_db')
  output_db=connect_db('neiva_output_db')
    
  if table_name=='processed ef':
    df=pd.read_sql(text('select * from Processed_EF'), con=output_db)
    efcoldf=pd.read_sql(text('select * from info_efcol_processed_data'), con=bk_db)

  if table_name=='integrated ef':
    df=pd.read_sql(text('select * from Integrated_EF'), con=output_db)
    efcoldf=pd.read_sql(text('select * from bkdb_info_efcol'), con=bk_db)

  iid=pcp.get_compounds(compound, 'name')[0].inchi
  ind=df[df['id']==iid].index[0]

  efcoldf[compound]=df[efcoldf['efcol']].iloc[ind].values

  ef_vals=list(efcoldf[compound][efcoldf['fire_type']==ft][efcoldf[compound].notna()])
  l1=efcoldf['study'][efcoldf['fire_type']==ft][efcoldf[compound].notna()]
  l2=efcoldf['fuel_type'][efcoldf['fire_type']==ft][efcoldf[compound].notna()]
  if table_name=='processed ef':
    ef_legend=list(l1)
  if table_name=='integrated ef':
    ef_legend=list(l1+':'+l2)

  ax = plt.subplot(111)
  plt.scatter(np.arange(len(ef_vals)), ef_vals)
  plt.ylabel('Emission factor (g/kg)')
  plt.title("Compound:"+compound+"; Fire type:"+ ft)
  plt.xticks(np.arange(len(ef_vals)), ef_legend, rotation=90)
  plt.grid(alpha=0.2)
  plt.tight_layout()
  return


def sum_ef_model_surrogate (chem, ft, model_surrogate):
    '''
    Calculates the VOC for a specific chemcial mechanism and fire type.
    '''
    output_db=connect_db('neiva_output_db')
    bk_db=connect_db('backend_db')
    
    # Loading Recommended_EF data, extracting NMOG_C and excluding compounds with 'unknown' (unknown chemical structure)
    dd=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
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
        ef=nmog_final['ef'][nmog_final[chem]==uu[i]]
        mm=nmog_final['mm'][nmog_final[chem]==uu[i]]
        
        totef=nmog_final['ef'][nmog_final[chem]==uu[i]].sum()
        weighted_mm=(ef*mm).sum()/ef.sum()
        prdf.loc[i,'ef']=totef
        prdf.loc[i,'weighted_mm']=weighted_mm
    
    # The TERP total EF is distributed equally to the following model species.
    if chem=='MOZT1':
        m_terp_sp=['MYRC','APIN','BPIN','LIMON']
        for i in range(len(m_terp_sp)):
            ind=prdf[prdf['MOZT1']==m_terp_sp[i]].index[0]
            prdf.loc[ind,'ef']=prdf['ef'].iloc[ind]+terp_unassigned/len(m_terp_sp)
            

    return prdf[[chem,'ef','weighted_mm']][prdf[chem]==model_surrogate]



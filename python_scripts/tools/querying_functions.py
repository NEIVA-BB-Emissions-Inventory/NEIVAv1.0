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

# This function returns the list of pollutant category of recommended ef dataset.
def display_pollutant_category():
    rdf=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
    ll=rdf['pollutant_category'].unique()
    return list(ll)

# This function returns the EF of specified pollutant category and fire type.
def select_ef_pollutant_category(ft, pc):
    efcol='AVG_'+ft.replace(' ','_')
    rdf=pd.read_sql(text('select * from Recommended_EF'), con=output_db)

    if pc=='PM optical property':
      return (rdf[['compound',efcol]][rdf['pollutant_category']==pc][rdf[efcol].notna()].reset_index(drop=True))
    else:
      return (rdf[['mm','formula','compound',efcol]][rdf['pollutant_category']==pc][rdf[efcol].notna()].reset_index(drop=True))

# This function returns the EF of specified compound name and table name (integrated ef, processed ef and recommended ef)
def select_compound(ft, com_name,table_name):
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

 # This funtion returns the model surrogates of a specified chemical mechanism
 # S07, S07T, S18B, S22, MOZT1
def model_surrogates(chem):
  pp=pd.read_sql(text('select * from Property_Surrogate'), con=output_db)
  return pp[chem].unique()

# This function returns the speciation compounds of specified chemical mechanism and model surrogate
def speciation_profile(ft,chem,spc):
  pp=pd.read_sql(text('select * from Property_Surrogate'), con=output_db)
  df=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
  df=df.merge(pp[pp.columns[3:]],on='id', how='left')
  efcol='AVG_'+ft.replace(' ','_')
  return df[['formula','smile',efcol,chem,'cstar','kOH']][df[chem]==spc][df[efcol].notna()]

# Plots ef data of a specified fire type and table name
def plot_ef(compound,ft, table_name):
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
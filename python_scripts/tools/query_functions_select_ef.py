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
from NEIVA.python_scripts.tools.assign_mozart_species import mozart_species
from NEIVA.python_scripts.data_integration_process.sort_molec_formula import *

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

    try:
        aa=pcp.get_compounds(com_name, 'name')[0].inchi
        ind=df[df['id']==aa].index[0]
        efcol=list(efcoldf['efcol'])

        efcoldf[com_name]=df[efcol].iloc[ind].values
        efcoldf[com_name]=round(efcoldf[com_name],3)
    
        ll=efcoldf[allcol][efcoldf['fire_type']==ft]
        ll=ll.sort_values(by='measurement_type')
        ll=ll[ll[com_name].notna()]
        ll=ll.reset_index(drop=True)
        return ll

    except:
        return 'Cannot assin ID. Search by formula'

def select_chemical_formula (ft, formula,table_name):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')

    if table_name=='integrated ef':
        df=pd.read_sql(text('select * from Integrated_EF'), con=output_db)
        efcoldf=pd.read_sql(text('select * from bkdb_info_efcol'), con=bk_db)

    if table_name=='processed ef':
        df=pd.read_sql(text('select * from Processed_EF'), con=output_db)
        efcoldf=pd.read_sql(text('select * from info_efcol_processed_data'), con=bk_db)
    
    ll=list(efcoldf['efcol'][efcoldf['fire_type']==ft])
    cols=['mm','formula','compound']+ll
    
    return df[cols][df['formula']==formula].reset_index(drop=True)



def select_compound_rdf (ft, com_name):
    output_db=connect_db('neiva_output_db')
    df=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
    
    try:
        aa=pcp.get_compounds(com_name, 'name')[0].inchi
        ind=df[df['id']==aa].index[0]
        col='AVG_'+ft.replace(' ','_')
        df[col]=round(df[col],3)
        return df[['mm','formula','compound',col]][ind:ind+1].reset_index(drop=True)
    except:
        return 'Cannot assign ID. Use chemical formula to search.'

def select_chemical_formula_rdf (ft, formula):
    output_db=connect_db('neiva_output_db')
    df=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
    
    col='AVG_'+ft.replace(' ','_')
    
    return df[['mm','formula', 'compound',col, 'id']][df['formula']==formula].reset_index(drop=True)
    


# Plots ef data of a specified fire type and table name
def get_ind (df, compound):
      if compound == 'PM<2.5':
          ind=list(df[df['pollutant_category']=='PM total'][df['compound'].str.contains('PM')].index)
          return ind
      if compound == 'BC':
          ind=list(df[df['id']=='BC'].index)
          return ind
      if compound == 'OC':
          ind=list(df[df['id']=='OC'].index)
          return ind
      if compound == 'NOx_as_NO':
          ind=list(df[df['id']=='NOx_as_NO'].index)
          return ind

      else:
          iid=pcp.get_compounds(compound, 'name')[0].inchi
          ind=list(df[df['id']==iid].index)
          return ind

def prepare_legend(fdf):
  for i in range(len(fdf)):
      fdf.loc[i,'legend']=fdf['study'].iloc[i]
      
  ii=list(fdf[fdf['fuel_type'].notna()].index)
  ii2=list(fdf[fdf['cookstove_name'].notna()].index)
  
  for i in ii:
      fdf.loc[i,'legend']=fdf['study'].iloc[i]+'\n'+fdf['fuel_type'].iloc[i]
  for i in ii2:
      fdf.loc[i,'legend']=fdf['study'].iloc[i]+'\n'+fdf['fuel_type'].iloc[i]+'\n'+fdf['cookstove_name'].iloc[i]
      
  return fdf

def plot_ef(compound,ft, table_name):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')
              
    if table_name=='processed ef':
        df=pd.read_sql(text('select * from Processed_EF'), con=output_db)
        efcoldf=pd.read_sql(text('select * from info_efcol_processed_data'), con=bk_db)
    
    if table_name=='integrated ef':
        df=pd.read_sql(text('select * from Integrated_EF'), con=output_db)
        efcoldf=pd.read_sql(text('select * from bkdb_info_efcol'), con=bk_db)

    try:
        iind=get_ind (df, compound)
        efcoldf[compound]=df[efcoldf['efcol']][df.index.isin(iind)].mean().values
        
        fdf=pd.DataFrame()
        fdf = efcoldf[efcoldf['fire_type']==ft][efcoldf[compound].notna()].reset_index(drop=True)
        fdf=prepare_legend(fdf)
        
        fdf=fdf.sort_values(by=compound, ascending=False)
        fdf=fdf.reset_index(drop=True)
        
        # Plot the figure   
        import seaborn as sns
        pal = sns.color_palette('bright',10)
        
        plt.figure(figsize=(12,8))
        ax1 = plt.subplot(111)
        
        x=np.arange(len(fdf))
        x_lab=list(fdf[compound][fdf['measurement_type']=='lab'].index)
        ef_lab=list(fdf[compound][fdf['measurement_type']=='lab'])
        
        plt.scatter(x, fdf[compound], zorder=3, color=pal[0], edgecolor='k', label='Field EF')
        plt.scatter(x_lab, ef_lab, zorder=3, color=pal[8], edgecolor='k', label='Lab EF')
        
        
        plt.ylabel('Emission factor (g/kg)', fontsize=11)
        
        plt.tick_params(labelsize=11)
        ax1.grid(linestyle='--',color='#EBE7E0',zorder=4)
        ax1.tick_params(axis='x',which='both',bottom=False)
        plt.setp(ax1.spines.values(),lw=1.5)
      
        plt.title("Compound: "+compound+"; Fire type:"+ ft, fontsize=11)
        plt.xticks(x, fdf['legend'], rotation=90)
        plt.legend(fontsize=10)
        plt.tight_layout()
    except:
         return 'Cannot assign ID. Use chemical formula to search.'
    return
     
def mce_vs_ef (compound, ft):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')
              
    df=pd.read_sql(text('select * from Integrated_EF'), con=output_db)
    efcoldf=pd.read_sql(text('select * from bkdb_info_efcol'), con=bk_db)

    try:
        iind=get_ind (df, compound)
        efcoldf[compound]=df[efcoldf['efcol']][df.index.isin(iind)].mean().values
        
        fdf=pd.DataFrame()
        if ft == 'all':
            fdf = efcoldf[efcoldf[compound].notna()].reset_index(drop=True)
        else: 
            fdf = efcoldf[efcoldf['fire_type']==ft][efcoldf[compound].notna()].reset_index(drop=True)
        fdf= fdf[fdf['MCE'].notna()]
        fdf=fdf.reset_index(drop=True)
        
        # Plot the figure   
        import seaborn as sns
        from sklearn.metrics import r2_score
        
        pal = sns.color_palette('bright',10)
        
        plt.figure(figsize=(10,6))
        ax1 = plt.subplot(111)
        
        x=np.arange(len(fdf))
        mce_lab=fdf['MCE'][fdf['measurement_type']=='lab']
        ef_lab=fdf[compound][fdf['measurement_type']=='lab']
        
        plt.scatter(fdf['MCE'], fdf[compound], zorder=3, color=pal[0], edgecolor='k', label='Field EF')
        plt.scatter(mce_lab, ef_lab, zorder=3, color=pal[8], edgecolor='k', label='Lab EF')
        
        plt.ylabel('Emission factor (g/kg)', fontsize=11)
        plt.xlabel('MCE', fontsize=11)
        
        if len(fdf)>4:
            coefficients = np.polyfit(fdf['MCE'], fdf[compound], 1)
            poly_fit = np.poly1d(coefficients)
            slope=round(coefficients[0],2)
            intercept=round(coefficients[1],2)
            # Calculate R-squared
            y_pred = poly_fit(fdf['MCE'])
            r_squared = r2_score(fdf[compound], y_pred)
            
            plt.plot(fdf['MCE'], poly_fit(fdf['MCE']), color=pal[6], label='y='+str(slope)+'*X + '+str(intercept) )
            #plt.scatter(0,0, color='k', label= 'R-squared:'+str(round(r_squared,2)) )            
        
        plt.tick_params(labelsize=11)
        ax1.grid(linestyle='--',color='#EBE7E0',zorder=4)
        ax1.tick_params(axis='x',which='both',bottom=False)
        plt.setp(ax1.spines.values(),lw=1.5)
      
        plt.title("Compound: "+compound+"; Fire type:"+ ft, fontsize=11)
        #plt.xticks(x, fdf['MCE'], rotation=90, fontsize=11)
        plt.legend(fontsize=11)
        plt.tight_layout()
    except:
         return 'Cannot assign ID. Use chemical formula to search.'
    return
         

def abundant_nmog (ft, chem, aa):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')
    
    efcol='AVG_'+ft.replace(' ','_')
    rdf=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
    rdf=rdf[rdf[efcol].notna()]
    rdf=rdf[rdf['pollutant_category']=='NMOC_g'].reset_index(drop=True)
    
    pp=pd.read_sql(text('select * from Property_Surrogate'), con=output_db)
    
    rdf=rdf.sort_values(by=efcol, ascending=False)
    rdf=rdf.reset_index(drop=True)
    
    rdf=rdf.merge(pp[['id',chem, aa]], on='id', how='left')
    
    return rdf[['mm','formula','compound', efcol, chem, aa, 'id']][:25].reset_index(drop=True)


def plot_abundant_nmog (ft):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')
    
    efcol='AVG_'+ft.replace(' ','_')
    rdf=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
    rdf=rdf[rdf[efcol].notna()]
    rdf=rdf[rdf['pollutant_category']=='NMOC_g'].reset_index(drop=True)
        
    rdf=rdf.sort_values(by=efcol, ascending=False)
    rdf=rdf.reset_index(drop=True)
    
    rdf=rdf[:25]
    
    import seaborn as sns
    pal = sns.color_palette('bright',10)
    
    plt.figure(figsize=(10,8))
    ax1 = plt.subplot(111)
    
    x=np.arange(len(rdf))
    
    plt.scatter(x, rdf[efcol], zorder=3, color=pal[0], edgecolor='k', label='Recommended EF')
    
    plt.ylabel('Emission factor (g/kg)', fontsize=11)
    plt.xlabel('Compound', fontsize=11)
            
    plt.tick_params(labelsize=11)
    ax1.grid(linestyle='--',color='#EBE7E0',zorder=4)
    ax1.tick_params(axis='x',which='both',bottom=False)
    plt.setp(ax1.spines.values(),lw=1.5)
  
    plt.title("Abundant NMOC_g' Fire type:"+ ft, fontsize=11)
    plt.xticks(x, rdf['compound'], rotation=90)
    plt.legend(fontsize=10)
    plt.tight_layout()
    
    return 

def nmog_with_high_n (ft, chem, aa):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')
    
    ncol='N_'+ft.replace(' ','_')
    efcol='AVG_'+ft.replace(' ','_')
    rdf=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
    rdf=rdf[rdf[efcol].notna()]
    rdf=rdf[rdf['pollutant_category']=='NMOC_g'].reset_index(drop=True)
    
    pp=pd.read_sql(text('select * from Property_Surrogate'), con=output_db)
    
    rdf=rdf.sort_values(by=ncol, ascending=False)
    rdf=rdf.reset_index(drop=True)
    
    rdf=rdf.merge(pp[['id',chem, aa]], on='id', how='left')
    
    return rdf[['mm','formula','compound', efcol, ncol,chem, aa, 'id']][:25].reset_index(drop=True)

def boxplot_abundant_nmog (ft):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')
    
    efcol='AVG_'+ft.replace(' ','_')
    rdf=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
    rdf=rdf[rdf[efcol].notna()]
    rdf=rdf[rdf['pollutant_category']=='NMOC_g'].reset_index(drop=True)
        
    rdf=rdf.sort_values(by=efcol, ascending=False)
    rdf=rdf.reset_index(drop=True)
    
    rdf=rdf[:25]
    
    pp=pd.read_sql(text('select * from Processed_EF'), con=output_db)
    efcoldf=pd.read_sql(text('select * from info_efcol_processed_data'), con=bk_db)
    efcol=efcoldf['efcol'][efcoldf['fire_type']==ft]
    
    import seaborn as sns
    pal = sns.color_palette('bright',10)
    
    plt.figure(figsize=(10,8))
    ax1 = plt.subplot(111)
    
    x=np.arange(len(rdf))
    
    for i in range(len(rdf)):
        aa=pp[efcol][pp['id']==rdf['id'].iloc[i]].T.dropna()
        vals=aa[aa.columns[0]].values   
        
        bp1 = ax1.boxplot(vals,showmeans=True,meanline=True,showfliers=True,patch_artist=True,positions=[i], widths=0.3,\
            medianprops=dict(linewidth=0),boxprops= dict(linewidth=1.5,edgecolor='k',facecolor=pal[0]),\
                  whiskerprops=dict(linestyle='-',linewidth=1,color='k'),\
                  meanprops=dict(color='red',linewidth=2,linestyle='-'),
                  flierprops = dict(marker='+',markerfacecolor=pal[8], markersize=7,))

    plt.yscale('log')
    plt.ylabel('Emission factor (g/kg)', fontsize=11)
    plt.xlabel('Compound', fontsize=11)
            
    plt.tick_params(labelsize=11)
    ax1.grid(linestyle='--',color='#EBE7E0',zorder=4)
    ax1.tick_params(axis='x',which='both',bottom=False)
    plt.setp(ax1.spines.values(),lw=1.5)
  
    plt.title("Abundant NMOC_g' Fire type:"+ ft, fontsize=11)
    plt.xticks(x, rdf['compound'], rotation=90)
    plt.legend(fontsize=10)
    plt.tight_layout()
    
    return 


def boxplot_ef (compound, ft, table_name):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')
              
    if table_name=='processed ef':
        df=pd.read_sql(text('select * from Processed_EF'), con=output_db)
        efcoldf=pd.read_sql(text('select * from info_efcol_processed_data'), con=bk_db)
    
    if table_name=='integrated ef':
        df=pd.read_sql(text('select * from Integrated_EF'), con=output_db)
        efcoldf=pd.read_sql(text('select * from bkdb_info_efcol'), con=bk_db)

    try:
        iind=get_ind (df, compound)
        efcol=efcoldf['efcol'][efcoldf['fire_type']==ft]
        aa=df[efcol][df.index.isin(iind)].T.dropna()
        vals=aa[aa.columns[0]].values   

        import seaborn as sns
        pal = sns.color_palette('bright',10)
        
        plt.figure(figsize=(5,8))
        ax1 = plt.subplot(111)
        x=[0]
        bp1 = ax1.boxplot(vals,showmeans=True,meanline=True,showfliers=True,patch_artist=True,positions=[i], widths=0.3,\
            medianprops=dict(linewidth=0),boxprops= dict(linewidth=1.5,edgecolor='k',facecolor=pal[0]),\
                  whiskerprops=dict(linestyle='-',linewidth=1,color='k'),\
                  meanprops=dict(color='red',linewidth=2,linestyle='-'),
                  flierprops = dict(marker='+',markerfacecolor=pal[8], markersize=7,))
    
        plt.ylabel('Emission factor (g/kg)', fontsize=11)
        plt.xlabel('Compound', fontsize=11)
                
        plt.tick_params(labelsize=11)
        ax1.grid(linestyle='--',color='#EBE7E0',zorder=4)
        ax1.tick_params(axis='x',which='both',bottom=False)
        plt.setp(ax1.spines.values(),lw=1.5)
      
        plt.title("Abundant NMOC_g' Fire type:"+ ft, fontsize=11)
        plt.xticks(x, compound, rotation=90)
        plt.legend(fontsize=10)
        plt.tight_layout()
    except:
         return 'Cannot assign ID. Use chemical formula to search.'
    
    return





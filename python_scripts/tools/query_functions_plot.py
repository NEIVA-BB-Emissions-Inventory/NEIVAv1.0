#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 10:48:36 2024

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
from NEIVA.python_scripts.tools.join_ef_property_table import *
from NEIVA.python_scripts.tools.query_functions_select_ef import get_ind

from sqlalchemy import text

def prepare_legend(fdf):
  for i in range(len(fdf)):
      fdf.loc[i,'legend']=fdf['study'].iloc[i]
  ii=list(fdf[fdf['fuel_type'].notna()].index)
  ii2=list(fdf[fdf['cookstove_name'].notna()].index)
  for i in ii:
      fdf.loc[i,'legend']=fdf['study'].iloc[i]+'\n'+'('+fdf['fuel_type'].iloc[i]+')'
  for i in ii2:
      fdf.loc[i,'legend']=fdf['study'].iloc[i]+'\n'+'('+fdf['fuel_type'].iloc[i]+'\n'+fdf['cookstove_name'].iloc[i]+')'
  return fdf

def plot_ef(ft, compound, table_name):
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
        
        plt.figure(figsize=(10,6))
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
         return 'Compound not found. Use chemical formula to search.'
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
            plt.text(0.8, 0.8, 'R-squared: '+str(round(r_squared,2)), fontsize=11, color='black', transform=plt.gca().transAxes)
        
        plt.tick_params(labelsize=11)
        ax1.grid(linestyle='--',color='#EBE7E0',zorder=4)
        ax1.tick_params(axis='x',which='both',bottom=False)
        plt.setp(ax1.spines.values(),lw=1.5)
        
        plt.text(0.8, 0.75, 'data count: '+str(len(fdf)), fontsize=11, color='black', transform=plt.gca().transAxes)
        #plt.text(0.8, 0.70, 'Compound: '+compound.capitalize(), fontsize=11, color='black', transform=plt.gca().transAxes)
        #plt.text(0.8, 0.65, 'Fire Type: '+ft.capitalize(), color='black', transform=plt.gca().transAxes)
  
        plt.title("Compound: "+compound.capitalize()+"; "+"Fire type: "+ ft.capitalize(), fontsize=11)
        #plt.xticks(x, fdf['MCE'], rotation=90, fontsize=11)
        plt.legend(fontsize=11, frameon=False)
        plt.tight_layout()
    except:
         return 'Compound not found. Use chemical formula to search.'
    return
         
def boxplot_abundant_nmog (ft):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')
    
    efcol='AVG_'+ft.replace(' ','_')
    ncol='N_'+ft.replace(' ','_')
    rdf=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
    rdf=rdf[rdf[efcol].notna()]
    rdf=rdf[rdf['pollutant_category']=='NMOC_g'].reset_index(drop=True)
        
    rdf=rdf.sort_values(by=efcol, ascending=False)
    rdf=rdf.reset_index(drop=True)
    
    rdf=rdf[:25]
    
    for i in range(len(rdf)):
        if len(rdf['compound'].iloc[i])>15:
            rdf.loc[i,'compound']=rdf['formula'].iloc[i]
    
    for i in range(len(rdf)):
        rdf.loc[i,'legend']=rdf['compound'].iloc[i]+';n='+str(rdf[ncol].iloc[i]).replace('.0','')

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
    #plt.xlabel('Compound', fontsize=11)
            
    plt.tick_params(labelsize=11)
    ax1.grid(linestyle='--',color='#EBE7E0',zorder=4)
    ax1.tick_params(axis='x',which='both',bottom=False)
    plt.setp(ax1.spines.values(),lw=1.5)
  
    plt.title("Abundant NMOC_g' Fire type:"+ ft, fontsize=11)
    plt.xticks(x, rdf['legend'], rotation=90)
    #plt.legend(fontsize=10)
    plt.tight_layout()
    
    return 

def boxplot_ef (compound, ft_list, table_name):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')
              
    if table_name=='processed ef':
        df=pd.read_sql(text('select * from Processed_EF'), con=output_db)
        efcoldf=pd.read_sql(text('select * from info_efcol_processed_data'), con=bk_db)
    
    if table_name=='integrated ef':
        df=pd.read_sql(text('select * from Integrated_EF'), con=output_db)
        efcoldf=pd.read_sql(text('select * from bkdb_info_efcol'), con=bk_db)
    
    if ft_list=='all':
        ft_list=['tropical forest','temperate forest','boreal forest',\
                 'savanna', 'crop residue', 'peat']
            
    import seaborn as sns
    pal = sns.color_palette('bright',10)
    x=np.arange(len(ft_list))
    plt.figure(figsize=(2.5*len(ft_list),7))
    ax1 = plt.subplot(111)
    
    legend=[]
    for i in range(len(ft_list)):
        iind=get_ind (df, compound)
        efcol=efcoldf['efcol'][efcoldf['fire_type']==ft_list[i]]
        aa=df[efcol][df.index.isin(iind)].mean().dropna().values
        vals=aa
        legend.append(ft_list[i].capitalize()+'(n='+str(len(aa))+')')
        
        bp1 = ax1.boxplot(vals,showmeans=True,meanline=True,showfliers=True,patch_artist=True,positions=[i], widths=0.3,\
            medianprops=dict(linewidth=0),boxprops= dict(linewidth=1.5,edgecolor='k',facecolor=pal[0]),\
                  whiskerprops=dict(linestyle='-',linewidth=1,color='k'),\
                  meanprops=dict(color='red',linewidth=2,linestyle='-'),
                  flierprops = dict(marker='+',markerfacecolor=pal[8], markersize=7,))
        
    plt.ylabel('Emission factor (g/kg)', fontsize=15)
    #plt.xlabel('Compound', fontsize=11)
            
    plt.tick_params(labelsize=15)
    ax1.grid(linestyle='--',color='#EBE7E0',zorder=4)
    ax1.tick_params(axis='x',which='both',bottom=False)
    plt.setp(ax1.spines.values(),lw=1.5)
      
    #plt.title("Compound: "+compound+"; "+"Fire type: "+ ft.capitalize(), fontsize=11)
    plt.xticks(x, legend, rotation=90, fontsize=15)
    #plt.legend(fontsize=10)
    plt.tight_layout()


def plot_model_surrogate (dd, ft, chem, model_surrogate):
    output_db=connect_db('neiva_output_db')
    bk_db=connect_db('backend_db')
    
    pp=pd.read_sql(text('select * from Processed_EF'), con=output_db)
    efcoldf=pd.read_sql(text('select * from info_efcol_processed_data'), con=bk_db)
    efcol=efcoldf['efcol'][efcoldf['fire_type']==ft]
    
    nmog=join_ef_property(dd)
    # Set EF column based on the input parameter 'fire type'
    avgcol='AVG_'+ft.replace(' ','_')
    ncol='N_'+ft.replace(' ','_')
    nmog['ef']=nmog[avgcol]
    nmog=nmog[nmog['ef'].notna()].reset_index(drop=True)
    
    nmog=nmog[nmog[chem]==model_surrogate]
    nmog=nmog.sort_values(by='ef', ascending=False)
    nmog=nmog.reset_index(drop=True)
    
    if len(nmog)>25:
        nmog=nmog[:25]
    
    for i in range(len(nmog)):
        if len(nmog['compound'].iloc[i])>15:
            nmog.loc[i,'compound']=nmog['formula'].iloc[i]
    
    for i in range(len(nmog)):
        legend=nmog['compound'].iloc[i]+';n='+str(nmog[ncol].iloc[i]).replace('.0','')
        nmog.loc[i,'legend']=legend

    import seaborn as sns
    pal = sns.color_palette('bright',10)
    
    plt.figure(figsize=(12,8))
    ax1 = plt.subplot(111)
    
    x=np.arange(len(nmog))
    
    for i in range(len(nmog)):
        aa=pp[efcol][pp['id']==nmog['id'].iloc[i]].T.dropna()
        vals=aa[aa.columns[0]].values   
        
        bp1 = ax1.boxplot(vals,showmeans=True,meanline=True,showfliers=True,patch_artist=True,positions=[i], widths=0.3,\
            medianprops=dict(linewidth=0),boxprops= dict(linewidth=1.5,edgecolor='k',facecolor=pal[0]),\
                  whiskerprops=dict(linestyle='-',linewidth=1,color='k'),\
                  meanprops=dict(color='red',linewidth=2,linestyle='-'),
                  flierprops = dict(marker='+',markerfacecolor=pal[8], markersize=7,))

    plt.yscale('log')
    plt.ylabel('Emission factor (g/kg)', fontsize=11)
    #plt.xlabel('Compound', fontsize=11)
    #plt.text(0.8, 0.9, 'Number of NMOC_g'+'('+model_surrogate+')'+': '+str(len(nmog_original)), fontsize=12, color='black', transform=plt.gca().transAxes)
    plt.tick_params(labelsize=11)
    ax1.grid(linestyle='--',color='#EBE7E0',zorder=4)
    ax1.tick_params(axis='x',which='both',bottom=False)
    plt.setp(ax1.spines.values(),lw=1.5)
  
    plt.title("Top NMOC_g with respect to EF", fontsize=12)
    plt.xticks(x, nmog['legend'], rotation=90)
    #plt.legend(fontsize=10)
    plt.tight_layout()
    
    return 
    



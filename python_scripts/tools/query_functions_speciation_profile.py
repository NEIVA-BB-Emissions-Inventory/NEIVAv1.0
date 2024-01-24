#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 20:02:31 2024

@author: samiha
"""

import pandas as pd
import numpy as np
import pubchempy as pcp
import matplotlib.pyplot as plt
from NEIVA.python_scripts.connect_with_mysql import*
from NEIVA.python_scripts.tools.assign_mozart_species import mozart_species
from NEIVA.python_scripts.data_integration_process.sort_molec_formula import *
from NEIVA.python_scripts.tools.join_ef_property_table import *

from sqlalchemy import text


# This function returns the speciation compounds of specified chemical mechanism and model surrogate
def speciation_profile(ft,chem,spc):
  bk_db=connect_db('backend_db')
  output_db=connect_db('neiva_output_db')
    
  pp=pd.read_sql(text('select * from Property_Surrogate'), con=output_db)
  pp=mozart_species(pp)
  
  df=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
  df=df.merge(pp[pp.columns[3:]],on='id', how='left')
  
  efcol='AVG_'+ft.replace(' ','_')
  return df[['mm','formula','compound','smile',efcol,chem]][df[chem]==spc][df[efcol].notna()]


def GFED_lumped_ef_calc (dd, chem, ft, model_surrogate):
    output_db=connect_db('neiva_output_db')
    bk_db=connect_db('backend_db')
    
    nmog=join_ef_property(dd)
    
    # Set EF column based on the input parameter 'fire type'
    efcol='AVG_'+ft.replace(' ','_')
    nmog['ef']=nmog[efcol]
    nmog=nmog[nmog['ef'].notna()].reset_index(drop=True)

    nmog_final = lump_com_with_speciation (nmog, chem)
    nmog_final = distribute_unk_ef (dd, efcol, nmog_final)
    
    nmog_final['ef']=nmog_final['ef']*nmog_final['conversion_factor']
    
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
    
    return prdf[[chem,'ef','weighted_mm']][prdf[chem]==model_surrogate]

def weighted_property (dd, ft, chem):
    output_db=connect_db('neiva_output_db')
    bk_db=connect_db('backend_db')
    
    nmog=join_ef_property(dd)
    # Set EF column based on the input parameter 'fire type'
    efcol='AVG_'+ft.replace(' ','_')
    nmog['ef']=nmog[efcol]
    nmog=nmog[nmog['ef'].notna()].reset_index(drop=True)

    nmog_final = lump_com_with_speciation (nmog, chem)
    nmog_final = distribute_unk_ef (dd, efcol, nmog_final)
    
    nmog_final['mole']=nmog_final['ef']/nmog_final['mm']
    # The list of unique model species
    uu=list(nmog_final[chem][nmog_final[chem].notna()].unique())
    
    # Creating the model species profile dataframe
    prdf=pd.DataFrame()
    for i in range(len(uu)):
        prdf.loc[i,chem]=uu[i]
        
        ef=nmog_final['ef'][nmog_final[chem]==uu[i]]
        mm=nmog_final['mm'][nmog_final[chem]==uu[i]]
        koh=nmog_final['kOH'][nmog_final[chem]==uu[i]]
        cstar=nmog_final['cstar'][nmog_final[chem]==uu[i]]
        hc=nmog_final['hc'][nmog_final[chem]==uu[i]]
        vp=nmog_final['vp'][nmog_final[chem]==uu[i]]
       
        totef=nmog_final['ef'][nmog_final[chem]==uu[i]].sum()
        totmole=nmog_final['mole'][nmog_final[chem]==uu[i]].sum()
        
        weighted_mm=(ef*mm).sum()/ef.sum()
        weighted_koh=(ef*koh).sum()/ef.sum()
        weighted_cstar=(ef*cstar).sum()/ef.sum()
        weighted_hc=(ef*hc).sum()/ef.sum()
        weighted_vp=(ef*vp).sum()/ef.sum()
      
        prdf.loc[i,'EF']=round(totef,3)
        prdf.loc[i,'Mole']=round(totmole,3)
        prdf.loc[i,'mm']=round(weighted_mm,3)
        prdf.loc[i,'kOH']=weighted_koh/1e-12
        prdf.loc[i,'cstar']=round(weighted_cstar,3)
        prdf.loc[i,'vp']=round(weighted_vp,3)
        prdf.loc[i,'hc']=weighted_hc
        
    prdf['Mole_fraction']=prdf['Mole']/prdf['Mole'].sum()
    prdf=prdf.sort_values(by='Mole_fraction', ascending=False)
    prdf=prdf.reset_index(drop=True)
     
    return prdf


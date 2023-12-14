#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 12:41:36 2022
@author: Samiha Shahid
"""

import pandas as pd
import numpy as np
import pubchempy as pcp
import sys

from utils import AltName,GrpCol

from connect_with_mysql import*
con=connect_db('NEIVA_db')

'''
Establishing Database Connections:
This section imports the necessary functions to connect to various 
databases and then initializes connections to 
five specific databases: NEIVA_db, legacy_db, raw_db, primary_db, and backend_db.
'''
from connect_with_mysql import *
n_con=connect_db('NEIVA_db')
legacy_db=connect_db('legacy_db')
raw_db=connect_db('raw_db')
primary_db=connect_db('primary_db')
bk_db=connect_db('backend_db')


def eliminate_general_terms(df):
    '''
    Compound name with general terms are eliminated. These terms typically refer to
    a class compound rather than specific compounds.
    Input:
        - DataFrame
    Return:
        -DataFrame: Filtered dataset excluding general compound names.
    '''
    g_terms=['isomer','other','various','general','DBE','C','substituted']
    for i in g_terms:
        df=df[~df['compound'].str.contains(i,na=False)]
    return df

def add_Spec2lumCom(lcdf,nmogdf):
    """
    Generates a speciation dataframe for lumped compounds.
    
    The function splits the lumped compounds of 'lcdf' at the '+' string and extracts individual compound names.
    Each name is assigned an InChI ID. The function then searches for these IDs in 'nmogdf' and appends them to 'lc_spec_df'.
    
    Input:
    - lcdf (DataFrame): Dataframe containing lumped compounds.
    - nmogdf (DataFrame): Dataframe with non-methane organic gases.
    
    Returns:
    - DataFrame: Speciation dataframe of lumped compounds.
    """    
    lc_spec_df=pd.DataFrame()
    for i in range(len(lcdf)):
        llid=[]
        df=pd.DataFrame()
        ll=lcdf['compound'][i].split('+')
        for k in ll:
            try:
                c=pcp.get_compounds(k.strip(),'name')
                llid.append(c[0].inchi)
                print(k.strip(), 'Assigned id')
            except:
                print(k.strip(), 'Unable to assign id.')
        specdf=nmogdf[nmogdf['id'].isin(llid)]
        if len(specdf)==len(ll):
            print('ids found in nmogdf')
            df=lcdf[i:i+1].append(specdf)
            lc_spec_df=lc_spec_df.append(df)
    lc_spec_df=lc_spec_df.reset_index(drop=True)
    return lc_spec_df


def Get_LumCom_Spec(nmogdf):
    """
    Extracts and aligns individual compounds from lumped compounds marked with a '+' sign.
    If found within the 'nmogdf', these individual compounds are then added to the 'lc_spec_df'.
    """
    # Identify compounds containing the '+' sign (indicating lumped compounds)
    com=[]
    for i  in range(len(nmogdf)):
        if nmogdf['compound'][i].find('+')!=-1:
            com.append(nmogdf['compound'].iloc[i])
    
    # Extract rows of lumped compounds from the dataset
    lcdf=nmogdf[nmogdf['compound'].isin(com)].reset_index(drop=True) 
    
    # Save to backend db for manual inspection
    lcdf[GrpCol(lcdf)[1]+['id']].to_sql(name='bkdb_nmog_LumpedCom',if_exists='replace',con=bk_db,index=False)
    
    # Get altered names to improve search accuracy in PubChem
    df_altName=pd.read_sql('select * from bkdb_nmog_LumpCom_altName',con=bk_db)
    lcdf=AltName(lcdf,df_altName)
    
   # Remove general terms from compound names
    lcdf=eliminate_general_terms(lcdf) 
    lcdf=lcdf.reset_index(drop=True)
    
    # Identify speciation for lumped compounds
    lc_spec_df=add_Spec2lumCom(lcdf,nmogdf) 
    
    return lc_spec_df





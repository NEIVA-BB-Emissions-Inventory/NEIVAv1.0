#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 12:41:36 2022
@author: Samiha Shahid
"""

import pandas as pd
import pubchempy as pcp
from sqlalchemy import text
from sqlalchemy import create_engine
import sys

from NEIVA.python_scripts.data_integration_process.data_formatting_functions import AltName,GrpCol
from NEIVA.python_scripts.data_integration_process.categorize_chemical_formula import *

'''
Establishing Database Connections:
This section imports the necessary functions to connect to various 
databases and then initializes connections to 
five specific databases: NEIVA_db, legacy_db, raw_db, primary_db, and backend_db.
'''
from NEIVA.python_scripts.connect_with_mysql import connect_db
# legacy_db=connect_db('legacy_db')
# raw_db=connect_db('raw_db')
# primary_db=connect_db('primary_db')
# bk_db=connect_db('backend_db')


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
    print('***************************************************************************************')
    print('This is a data sorting step. The following steps are executed:')
    print('1. Split a lumped compound into individual components, assign id to the componenets')
    print('2. Search the ids within the integrated dataset')
    print('3. Align the lumped compound and individual compounds if found.')
    print('***************************************************************************************')
 
    lc_spec_df=pd.DataFrame()
    for i in range(len(lcdf)):
        llid=[]
        df=pd.DataFrame()
        ll=lcdf['compound'][i].split('+')
        print('Lumped compound- '+lcdf['compound'][i])
       
        for k in ll:
            try:
                c=pcp.get_compounds(k.strip(),'name')
                llid.append(c[0].inchi)
                print('Assigned id to individual compound: ', k.strip())
            except:
                print('Unable to assign id to inidividual compound: ', k.strip())
        specdf=nmogdf[nmogdf['id'].isin(llid)]
        if len(specdf)==len(ll):
            print('All individual ids are found in the Integrated Dataset')
            # df=lcdf[i:i+1].append(specdf)
            df = pd.concat([lcdf[i:i+1], specdf], ignore_index=True)
            # lc_spec_df=lc_spec_df.append(df)
            lc_spec_df = pd.concat([lc_spec_df,df], ignore_index=True)
        print('__________________________________________________________________')
    lc_spec_df=lc_spec_df.reset_index(drop=True)
    return lc_spec_df


def sync_lumped_compound_and_speciation(nmogdf):
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
    engine = create_engine("mysql+pymysql://root:root@localhost/"+'backend_db')
    lcdf[GrpCol(lcdf)[1]+['id']].to_sql(name='bkdb_nmog_LumpedCom',if_exists='replace',con=engine,index=False)
    
    # Get altered names to improve search accuracy in PubChem
    bk_db=connect_db('backend_db')
    df_altName=pd.read_sql(text('select * from bkdb_nmog_LumpCom_altName'),con=bk_db)
    lcdf=AltName(lcdf,df_altName)
    
   # Remove general terms from compound names
    lcdf=eliminate_general_terms(lcdf) 
    lcdf=lcdf.reset_index(drop=True)
    
    # Identify speciation for lumped compounds
    lc_spec_df=add_Spec2lumCom(lcdf,nmogdf) 
    
    return lc_spec_df

def import_fc_dataset(nmogdf,lc_spec_df):
    '''
    Imports specific and simple fractional contribution datasets to the backend database.
    
    The function processes two dataframes, nmogdf and lc_spec_df, to derive two resulting datasets: 
    specific_fc_df and simple_fc. The specific dataset contains aligned lumped compounds and speciation,
    while the simple dataset consists of single lumped compounds for each formula.

    Parameters:
    - nmogdf: Dataframe with total NMOG data.
    - lc_spec_df: Dataframe with lumped compounds and their speciation.
    
    Returns:
    None. The function performs in-place modifications and exports results to the database.
    '''    
    # Construct specific fractional contribution dataset with aligned lumped compounds and speciation.
    engine = create_engine("mysql+pymysql://root:root@localhost/"+'backend_db')
    specific_fc_df=lc_spec_df[GrpCol(lc_spec_df)[1]+['id','study']]
    specific_fc_df.to_sql(name='bkdb_fc_calc_specific', con=engine, if_exists='replace', index=False)
    
    # Extract simple fractional contribution dataset with only one lumped compound per formula.
    simple_fc=nmogdf[nmogdf['formula'].isin(assign_formula_type(nmogdf)[3])]
    # Exclude entries present in specific_fc_df.
    simple_fc=simple_fc[~simple_fc['formula'].isin(specific_fc_df['formula'].tolist())]
    #____Excluding Isoprene and Toluene______________________________________________
    iid=['InChI=1S/C5H8/c1-4-5(2)3/h4H,1-2H2,3H3',
         'InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3']
    simple_fc=simple_fc[~simple_fc['id'].isin(iid)]
    # Rearrange columns and export to backend database.
    simple_fc=simple_fc[GrpCol(simple_fc)[1]+['id','study']]
    simple_fc.to_sql(name='bkdb_fc_calc_simple',con=engine, if_exists='replace', index=False)
    print('Imported fractional contribution datasets to Backend database-')
    return



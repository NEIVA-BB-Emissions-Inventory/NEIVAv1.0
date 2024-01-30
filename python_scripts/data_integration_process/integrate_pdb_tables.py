"""
Created on Tue Mar  8 12:07:36 2022
@author: Samiha Shahid
"""

import pandas as pd
from sqlalchemy import text
from NEIVA.python_scripts.data_integration_process.sort_molec_formula import *
from NEIVA.python_scripts.data_integration_process.display_pretty_table import *

'''
Establishing Database Connections:
This section imports the necessary functions to connect to various 
databases and then initializes connections to 
five specific databases: NEIVA_db, legacy_db, raw_db, primary_db, and backend_db.
'''
from NEIVA.python_scripts.connect_with_mysql import connect_db, get_table_name
# legacy_db=connect_db('legacy_db')
# raw_db=connect_db('raw_db')
# primary_db=connect_db('primary_db')
# bk_db=connect_db('backend_db')


def integrate_tables():
    '''
    Integrates all primary_db datasets into a single dataframe.
    
    - Fetches data from primary database tables.
    - Drops unnecessary columns and ensures dataset consistency.
    - Asserts the absence of duplicate IDs and checks for rows with all NaN values.
    - Processes the formula and calculates the nominal molar mass.
    
    Returns:
        df: Integrated dataframe.
    '''    
    df=pd.DataFrame()
    idcols=['mm','formula','compound','pollutant_category','id']
    
    #__ table names of primary database
    primary_db = connect_db('primary_db')
    pdb_tbl_names = get_table_name('primary_db')
    
    df=pd.read_sql(text('select * from '+'pdb_koss18'),con=primary_db)
    dropcol = df.columns[~df.columns.isin(idcols+df.filter(like='EF').columns.tolist()+['id'])].tolist()
    df=df.drop(columns=dropcol)    
    
    tbl_names= set(pdb_tbl_names)-{'pdb_koss18'}
    
    tbldf=pd.DataFrame()
    tbldf['tbl_names']=list(tbl_names)
    
    tbldf=tbldf.sort_values(by='tbl_names')
    tbldf=tbldf.reset_index(drop=True)
    
    print(tbldf['tbl_names'])
    
    for i in range(len(tbldf)):
        data=pd.read_sql(text('select * from '+tbldf['tbl_names'].iloc[i]),con=primary_db)
        data_id=data[idcols]
        data_ef=data[data.filter(like='EF').columns.tolist()+['id']]
        unmatched=data_id[~data_id['id'].isin(df['id'])]
        # df=df.append(unmatched)
        df = pd.concat([df,unmatched], ignore_index=True)
        df = df.merge(data_ef,on='id',how='left')
        
    assert len(df[df['id'].duplicated()]) == 0
    assert len(df)==len(df[df.columns[df.columns.str.contains('EF',na=False)]].dropna(how='all'))
    
    df=order_formula(df)
    df=nominal_mm_calulator(df)
    
    t=create_PrettyTable_col2(['Check','Output'],['Duplicate IDs','Row with all NaNs','formula column'],['Clear','Clear','Cleaned'])
    print('\nData Integration process COMPLETE!')
    print(t)
    print('Length of Integrated dataset :'+str(len(df)))
    return df


def sort_inorganic_gas_data(df):
    '''
    Sorts the input dataset by molar mass for inorganic gases and methane.
    Input:
        df (DataFrame): Input dataset with a 'pollutant_category' column.
    Returns:
        DataFrame: Sorted subset of the input dataset.
    '''    
    igdf=df[df['pollutant_category'].isin(['inorganic gas','methane'])]
    igdf=igdf.sort_values(by='mm').reset_index(drop=True)
    return igdf

def sort_particulate_matter_data(df):
    '''
    Sorts the particulate matter dataset based on a predefined sequence.
    Input:
        df (DataFrame): Input dataset with a 'pollutant_category' column.
    Returns:
        DataFrame: Ordered subset of the input dataset.
    Note:
        The order sequence is retrieved from the 'bkdb_pm_order_seq' table in 'backend_db'.
    '''    
    bk_db=connect_db('backend_db')
    pm_arrange_seq=pd.read_sql(text('select * from bkdb_pm_order_seq'),con=bk_db)
    pm_arrange_seq=list(pm_arrange_seq['pollutant_category_p'])
    
    pmdf=pd.DataFrame()
    for pm_type in pm_arrange_seq:
        # pmdf=pmdf.append(df[df['pollutant_category']==i])
        pmdf = pd.concat([pmdf,df[df['pollutant_category']==pm_type]], ignore_index=True)
    
    return pmdf

def select_nmog(df):
    nmogdf=df[df['pollutant_category']=='NMOC_g']
    nmogdf=nmogdf.sort_values(by=['mm','formula','id']).reset_index(drop=True)
    print('Non-Methane Organic Compounds Gas-Phase (NMOC_g) Data Frame: [ROW, COLUMN] ='+ '['+str(len(nmogdf))+' '+str(len(nmogdf.columns))+']' )
    return nmogdf

def sort_nmog_data(nmogdf):
    '''
    Sorts the NMOG (Non-Methane Organic Gases) dataframe based on molecular weight and 
    different criteria of the 'id' column.
    
    The function primarily sorts the dataframe first by the molecular mass. Within each molecular mass, 
    it further sorts rows based on the presence of 'InChI' in the 'id' column, followed by 
    rows without 'InChI' but not in the rdb_hatch15 dataset, and then rows with ids present 
    in the rdb_hatch15 dataset.
    
    Parameters:
    - nmogdf: Input dataframe containing NMOG data.
    
    Returns:
    - nmogdf_sorted: A sorted dataframe based on the above criteria.
    '''
    # Loadig the 'rdb_hatch15' dataste.
    raw_db=connect_db('raw_db')
    hid = pd.read_sql(text('select * from rdb_hatch15'),con=raw_db)
    hid = hid[hid['h_id'].notna()].reset_index(drop=True)
    
    # Sorting nmogdf by molecular mass.
    nmogdf=nmogdf.sort_values(by=['mm']).reset_index(drop=True)
    
    # List of unique formula of 'nmogdf'
    uf=list(nmogdf['formula'].unique())
    
    # Sorting each unique formula based on the criteria.
    nmogdf_sorted=pd.DataFrame()
    for f in uf:
        formula_df = nmogdf[nmogdf['formula']==f]
        # Select rows where 'id' contains 'InChI'
        aa1=formula_df[formula_df['id'].str.contains('InChI')]
        aa1=aa1.sort_values(by='id')
        
        # Select rows where 'id' does not contains InChI and is not in hid.
        aa2=formula_df[~formula_df['id'].str.contains('InChI')][~formula_df['id'].isin(hid['h_id'])]
        
        # Select rows where 'id' does not contains InChI but is in hid.
        aa3=formula_df[~formula_df['id'].str.contains('InChI')][formula_df['id'].isin(hid['h_id'])]
        
        # aa=aa2.append(aa1).append(aa3)
        aa = pd.concat([aa2,aa1,aa3], ignore_index=True)
        # nmogdf_sorted=nmogdf_sorted.append(aa).reset_index(drop=True)
        nmogdf_sorted = pd.concat([nmogdf_sorted,aa], ignore_index=True)
    
    return nmogdf_sorted

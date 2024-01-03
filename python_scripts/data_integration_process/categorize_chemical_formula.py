#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  1 13:36:48 2024

@author: samiha
"""

import pandas as pd
from sqlalchemy import text
from NEIVA.python_scripts.connect_with_mysql import connect_db

# primary_db=connect_db('primary_db')

def assign_formula_type(df):
    '''
    Categorizes chemical formulas from the input dataframe.
    
    The categorization is done based on formula occurrence and the presence of lumped compounds.
    A lumped compound is determined if its ID doesn't match the pdb_hatch15 table or an InChI structure.

    Parameters:
    - df: Dataframe containing chemical data.

    Returns:
    - f_no_spec: Formulas without speciation.
    - f_spec_no_lc: Formulas with speciation but without lumped compounds.
    - f_spec_multiple_lc: Formulas with speciation and multiple lumped compounds.
    - f_spec_one_lc: Formulas with speciation and a single lumped compound.
    
    '''
    # Loading the hatch15 from primary database.
    primary_db=connect_db('primary_db')
    hid = pd.read_sql(text('select * from pdb_hatch15'),con=primary_db)
    
    # Unique formula of 'df' dataframe.
    all_formula=df['formula'].unique().tolist()
    
    # Identify formula that only appear once in the dataframe.
    f_no_spec=[]
    for f in all_formula:
        if len(df[df['formula']==f])==1: 
            f_no_spec.append(f)
    
    # Formula with more than one occurance.
    f_spec=set(all_formula)-set(f_no_spec)
    
    # Formulas with speciation and without lumped compounds.
    f_spec_no_lc=[]
    for f in f_spec:
        aa=df[df['formula']==f]
        bb=aa[aa['id'].str.contains('InChI',na=False)]
        if len(aa)==len(bb):
            f_spec_no_lc.append(f)
    
    # Formulas that are not in 'f_spec_no_lc'
    f_spec_lc=f_spec-set(f_spec_no_lc)
    
    # Formulas with multiple lumped compound.
    f_spec_multiple_lc=[]
    for f in f_spec_lc:
        aa=df[df['formula']==f]
        aa=aa[~aa['id'].str.contains('InChI',na=False)]
        aa=aa[~aa['id'].isin(hid['id'].tolist())]
        if len(aa)>1:
             f_spec_multiple_lc.append(f)
    # Formulas with a single lumped compound.
    f_spec_one_lc=[]
    for f in f_spec_lc:
        aa=df[df['formula']==f]
        aa=aa[~aa['id'].str.contains('InChI',na=False)]
        aa=aa[~aa['id'].isin(hid['id'])]
        if len(aa)==1:
            f_spec_one_lc.append(f)
    
    return f_no_spec, f_spec, f_spec_multiple_lc, f_spec_one_lc 

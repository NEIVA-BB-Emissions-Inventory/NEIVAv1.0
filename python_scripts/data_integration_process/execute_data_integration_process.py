#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Samiha Binte Shahid
Email: sbint003@ucr.edu
Github: @SamihaShahid

"""
import pandas as pd
import numpy as np

from NEIVA.python_scripts.data_integration_process.integrate_pdb_tables import *
from NEIVA.python_scripts.data_integration_process.merge_lumped_com import *
from NEIVA.python_scripts.data_integration_process.align_lumped_com_and_spec import *
from NEIVA.python_scripts.data_integration_process.data_formatting_functions import  *
from NEIVA.python_scripts.connect_with_mysql import connect_db, get_table_name

output_db=connect_db('neiva_output_db')


'''
This script integrates primary database tables into a unified dataset. 
The functions employed in this integration
process are sourced from DataInt_mainFunc.py, Merge_MultLumpCom.py,
Get_LumpCom_Spec.py, and utils.py scripts. 
The final, integrated dataset is saved as IntData in the output_db.
'''

# Integrates the primary database datasets.
int_df=integrate_tables()

# Extracts dataframe where pollutant_category is NMOC_g 
nmogdf=select_nmog(int_df)   
nmogdf=assign_study_column(nmogdf)

# Processes and manages record with multiple lumped compounds.
r_iddf, iddf =merge_lumped_compound_same_formula(nmogdf)

# Updates nmogdf by replacing iddf with r_iddf
nmogdf=insert_rdf_nmogdf(nmogdf,r_iddf,iddf) 

# Decomposes lumped compounds and aligns with individual compounds
lc_spec_df=sync_lumped_compound_and_speciation(nmogdf)

# Sorts nmogdf for further processing
nmogdf=sort_nmog_data(nmogdf)

# Loads the dataset for fractional contribution calculations to the backend database
import_fc_dataset(nmogdf,lc_spec_df)

# Sorts and merges inorganic gases, particulate matter data with nmogdf
igdf=sort_inorganic_gas_data(int_df)
pmdf=sort_particulate_matter_data(int_df)
mdf = pd.concat([igdf, nmogdf, pmdf], ignore_index=True)
# Storing the final integrated dataset in the database 
mdf=mdf.reset_index(drop=True)

# Adjusts datatype for the 'mm' column
mdf=str_float(mdf,'mm')

mdf.to_sql(name='Integrated_EF',con=output_db, if_exists='replace', index=False)
print('Updated the Integrated dataset in NEIAV_db.sql')


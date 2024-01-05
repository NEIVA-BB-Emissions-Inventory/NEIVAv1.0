#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 19:24:21 2022

@author: Samiha Shahid
"""
import pandas as pd

from NEIVA.python_scripts.data_processing_steps.data_calculations import *
from NEIVA.python_scripts.data_processing_steps.lab_data_emission_ratio_adjust import * 
from NEIVA.python_scripts.data_processing_steps.assign_fractional_contribution import *
from NEIVA.python_scripts.data_integration_process.data_formatting_functions import  GrpCol, rearrange_col_finaldf

from NEIVA.python_scripts.data_processing_steps.info_table_sorting_functions import assign_year_col_efcoldf, assign_legend_col

# Connect to MySQL database
from NEIVA.python_scripts.connect_with_mysql import connect_db
output_db=connect_db('neiva_output_db')
bk_db=connect_db('backend_db')

# Read EF column information from 'bkdb_info_efcol' table.
efcoldf=pd.read_sql('select * from bkdb_info_efcol', con=bk_db) 

# Load the integrated datase.
intdf=pd.read_sql('select * from Integrated_EF', con=output_db) 

# Calculate NOx as equivalent to NO.
intdf=calc_NOx_as_NO (intdf)

# Calculate lab study and update EF column information table.
intdf_2, efcoldf = calculate_average_lab_study(intdf, efcoldf)

# Check the statement to ensure data consistency.
assert len(intdf_2.columns[intdf_2.columns.str.contains('EF')])== len(efcoldf)

# Assign year and legend columns in EF column information table.
efcoldf = assign_year_col_efcoldf(efcoldf)
efcoldf = assign_legend_col(efcoldf)

# Perform emission ratio adjustment calculations on the integrated dataset.
intdf_3=lab_data_adjust_to_field_conditions(intdf_2,efcoldf)[0]

# Assign 'N' columns to the integrated dataset based on EF column information.
intdf_3=assign_data_count_column(intdf_3,efcoldf)

# Calculate average values for the integrated dataset.
avgdf = calculate_average_fire_types(intdf_3, efcoldf)
# Round average columns to four decimal places.
avgdf = round_avg_cols(avgdf)

# Calculate fractional contributions in the average dataset.
avgdf=calculalate_fractional_contribution(avgdf)  

# Rearrange columns in the final dataset.
#avgdf=rearrange_col_finaldf(avgdf)

# Store the Recommended EF dataset in the 'Recommended_EF' table of 'neiva_output_db'.
avgdf.to_sql(name = 'Recommended_EF', con=output_db, if_exists='replace', index=False)


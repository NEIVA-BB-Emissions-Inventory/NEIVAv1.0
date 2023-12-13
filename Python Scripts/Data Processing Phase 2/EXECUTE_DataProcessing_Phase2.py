#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 19:24:21 2022

@author: Samiha Shahid
"""
import pandas as pd

from utils_calc import calc_NOx_as_NO
from ER_ADJ_calc import * 
from AVG_n_FC_calc import assign_n_cols,get_avg_df,fc_calc, round_avg_cols
from utils import GrpCol, rearrange_col_finaldf
from utils_sort import assign_year_col_efcoldf, assign_legend_col

# Connect to MySQL database
from connect_with_mysql import connect_db
output_db=connect_db('neiva_output_db')
bk_db=connect_db('backend_db')

# Read EF column information from 'bkdb_info_efcol' table.
efcoldf=pd.read_sql('select * from bkdb_info_efcol', con=bk_db) 

# Load the integrated datase.
intdf=pd.read_sql('select * from IntData', con=output_db) 

# Calculate NOx as equivalent to NO.
intdf=calc_NOx_as_NO (intdf)

# Calculate lab study and update EF column information table.
intdf_2, efcoldf = get_lab_study_avg(intdf, efcoldf)

# Check the statement to ensure data consistency.
assert len(intdf_2.columns[intdf_2.columns.str.contains('EF')])== len(efcoldf)

# Assign year and legend columns in EF column information table.
efcoldf = assign_year_col_efcoldf(efcoldf)
efcoldf = assign_legend_col(efcoldf)

# Perform emission ratio adjustment calculations on the integrated dataset.
intdf_3=er_adj(intdf_2,efcoldf)[0]

# Assign 'N' columns to the integrated dataset based on EF column information.
intdf_3=assign_n_cols(intdf_3,efcoldf)

# Calculate average values for the integrated dataset.
avgdf = get_avg_df(intdf_3, efcoldf)
# Round average columns to four decimal places.
avgdf = round_avg_cols(avgdf)

# Calculate fractional contributions in the average dataset.
avgdf=fc_calc(avgdf)  

# Rearrange columns in the final dataset.
avgdf=rearrange_col_finaldf(avgdf)

# Store the Recommended EF dataset in the 'Recommended_EF' table of 'neiva_output_db'.
avgdf.to_sql(name = 'Recommended_EF', con=output_db, if_exists='replace', index=False)


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 13:21:31 2025

@author: samiha
"""
from utils_calc import calc_NOx_as_NO
from ER_ADJ_calc import * 
from AVG_n_FC_calc import assign_n_cols,get_avg_df,fc_calc, round_avg_cols, exception_cases
from utils import GrpCol, rearrange_col_finaldf
from utils_sort import assign_year_col_efcoldf, assign_legend_col

from down_weight_travis_savanna import down_weight_travis_savanna 

# Connect to MySQL database
from connect_with_mysql import connect_db
output_db=connect_db('neiva_output_db')
bk_db=connect_db('backend_db')
primary_db=connect_db('primary_db')


pdb_tbl_names = get_table_name('primary_db')

# ________________________Fix BC of boreal forest_____________________________
df=pd.read_sql('select * from '+'pdb_bf_hayden22',con=primary_db)
bc_ef = df['EF_boreal_forest_hayden22'][df['id']=='BC'].values[0]
pm1_ef = df['EF_boreal_forest_hayden22'][df['id']=='PM1'].values[0]

r=bc_ef/pm1_ef

#_____________________________________________________________________________

# Read EF column information from 'bkdb_info_efcol' table.
efcoldf=pd.read_sql('select * from bkdb_info_efcol', con=bk_db) 
# Load the integrated datase.
intdf=pd.read_sql('select * from Integrated_EF', con=output_db) 

# __________ fix savanna _____________________________________________________
# drop EF_savanna_desservettaz20 and EF_grassland_travis23
intdf = intdf.drop(['EF_savanna_desservettaz20','EF_grassland_travis23'], axis=1)
efcoldf = efcoldf[~efcoldf['efcol'].isin(['EF_savanna_desservettaz20', 'EF_grassland_travis23'])]
efcoldf=efcoldf.reset_index(drop=True)
#____________________________________________________________________________-
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

#________________ estimate new BC EF boreal forest ___________________________
bf_ef = efcoldf['efcol'][efcoldf['fire_type']=='boreal forest']
pm_2_5_all = intdf_3[bf_ef][intdf_3['pollutant_category']=='PM total'].mean().values

pm_2_5_all_corrected = pm_2_5_all * r

BC_EF = np.nanmean(pm_2_5_all_corrected)

#____________________________________________________________________-
# Calculate average values for the integrated dataset.
avgdf = get_avg_df(intdf_3, efcoldf)
# Round average columns to four decimal places.
#avgdf = round_avg_cols(avgdf)

# Calculate fractional contributions in the average dataset.
avgdf=fc_calc(avgdf)  

# Rearrange columns in the final dataset.
#avgdf=rearrange_col_finaldf(avgdf)

# Store the Recommended EF dataset in the 'Recommended_EF' table of 'neiva_output_db'.
# add the exceltion cases
avgdf=exception_cases(avgdf)

# update in avgdf
bc_ind=avgdf[['compound','AVG_boreal_forest']][avgdf['compound']=='BC'].index[0]
avgdf.loc[bc_ind, 'AVG_boreal_forest']=BC_EF

# avgdf.to_excel('/Users/samiha/Desktop/NEIVA_final_v2/Recommended_EF_GFED.xlsx',index=False)
avgdf.to_sql(name = 'Recommended_EF', con=output_db, if_exists='replace', index=False)



# avgdf.to_csv('/Users/samiha/Desktop/git_package/NEIVA/data/CSV/Recommended_EF.csv', index=False)












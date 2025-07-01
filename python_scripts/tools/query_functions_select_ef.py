#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 12:34:51 2024

@author: samiha
"""

import pandas as pd
import numpy as np
import pubchempy as pcp
import matplotlib.pyplot as plt
from NEIVA.python_scripts.connect_with_mysql import connect_db
from NEIVA.python_scripts.tools.assign_mozart_species import mozart_species
from NEIVA.python_scripts.data_integration_process.sort_molec_formula import *
from NEIVA.python_scripts.tools.assign_geos_chem_species import geos_chem_species
from NEIVA.python_scripts.tools.number_format_function import rounding
from NEIVA.python_scripts.tools.join_ef_property_table import *

#from NEIVA.python_scripts.tools.query_functions_plot import get_ind


from sqlalchemy import text


def get_ind(df, compound):
    if compound == 'PM2.5*':
        ind = df[df['compound'].str.contains('PM')][df['id']!='PM10'][df['id']!='PM2.5_ipcc'].index.tolist()
    elif compound in ['PM10', 'PM2.5','PM1','PM2.5(PM1-PM5)','OA', 'EC', 'BC', 'OC', 'NOx_as_NO']:
        ind = df[df['id'] == compound].index.tolist()
    else:
        iid = pcp.get_compounds(compound, 'name')[0].inchi
        ind = df[df['id'] == iid].index.tolist()
    return ind

def get_ind_rdf(df, compound):
    if compound in ['PM2.5*', 'PM10', 'OA', 'EC', 'BC', 'OC', 'NOx_as_NO']:
        ind = df[df['id'] == compound].index.tolist()
    else:
        iid = pcp.get_compounds(compound, 'name')[0].inchi
        ind = df[df['id'] == iid].index.tolist()
    
    return ind

# This function returns the PM2.5, OC, BC data for a specified fire type and table name. The table name inlcudes three different tables
# integrated ef, processed ef, recommended ef.
def select_pm_data (ft, table_name):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')
    if table_name=='integrated ef':
        efcoldf=pd.read_sql(text('select * from bkdb_info_efcol'), con=bk_db)
        df=pd.read_sql(text('select * from Integrated_EF'), con=output_db)
        
        efcol=list(efcoldf['efcol'])
        efcoldf['BC']=df[efcol][df['id']=='BC'].values[0]
        efcoldf['OC']=df[efcol][df['id']=='OC'].values[0]
        efcoldf['PM2.5']=df[efcol][df['id']=='PM2.5'].values[0]
        efcoldf['PM10']=df[efcol][df['id']=='PM10'].values[0]
        efcoldf['PM1']=df[efcol][df['id']=='PM1'].values[0]
        efcoldf['PM2.5(non-pile burning)']=df[efcol][df['id']=='PM2.5_npb'].values[0]
        efcoldf['PM2.5(pile-burning)']=df[efcol][df['id']=='PM2.5_pb'].values[0]
        efcoldf['OA']=df[efcol][df['id']=='OA'].values[0]
        efcoldf['EC']=df[efcol][df['id']=='EC'].values[0]
        if ft=='crop residue':
            fdf=efcoldf[['legend','MCE','PM2.5','PM1','PM10','PM2.5(non-pile burning)','PM2.5(pile-burning)','OC','BC','OA', 'EC']][efcoldf['fire_type']==ft].reset_index(drop=True)
            fdf=fdf.applymap(lambda x: rounding(x))
            return fdf
        if ft!='crop residue':
            fdf=efcoldf[['legend','MCE','PM2.5','PM1','PM10','OA','OC','BC','EC']][efcoldf['fire_type']==ft].reset_index(drop=True)
            fdf=fdf.applymap(lambda x: rounding(x))
            return fdf
    if table_name=='processed ef':
        efcoldf=pd.read_sql(text('select * from info_efcol_processed_data'), con=bk_db)
        df=pd.read_sql(text('select * from Processed_EF'), con=output_db)
        
        efcol=list(efcoldf['efcol'])
        efcoldf['BC']=df[efcol][df['id']=='BC'].values[0]
        efcoldf['OC']=df[efcol][df['id']=='OC'].values[0]
        efcoldf['PM2.5']=df[efcol][df['id']=='PM2.5'].values[0]
        efcoldf['PM10']=df[efcol][df['id']=='PM10'].values[0]
        efcoldf['PM1']=df[efcol][df['id']=='PM1'].values[0]
        efcoldf['PM2.5(non-pile burning)']=df[efcol][df['id']=='PM2.5_npb'].values[0]
        efcoldf['PM2.5(pile-burning)']=df[efcol][df['id']=='PM2.5_pb'].values[0]
        efcoldf['OA']=df[efcol][df['id']=='OA'].values[0]
        efcoldf['EC']=df[efcol][df['id']=='EC'].values[0]
        if ft=='crop residue':
            fdf=efcoldf[['legend','MCE','PM2.5','PM1','PM10','PM2.5(non-pile burning)','PM2.5(pile-burning)','OC','BC','OA', 'EC']][efcoldf['fire_type']==ft].reset_index(drop=True)
            fdf=fdf.applymap(lambda x: rounding(x))
            return fdf
        if ft!='crop residue':
            fdf=efcoldf[['legend','MCE','PM2.5','PM1','PM10','OA','OC','BC','EC']][efcoldf['fire_type']==ft].reset_index(drop=True)
            fdf=fdf.applymap(lambda x: rounding(x))
            return fdf
    if table_name=='recommended ef':
        df=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
        efcol='AVG_'+ft.replace(' ','_')
        iid=['PM2.5*','PM10','OA','OC','BC','EC']
        fdf=df[['compound',efcol]][df['id'].isin(iid)].reset_index(drop=True)
        fdf=fdf.applymap(lambda x: rounding(x))
        return fdf


# This function returns the EF of specified pollutant category and fire type.
def select_ef_pollutant_category(ft, pc):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')
    efcol='AVG_'+ft.replace(' ','_')
    rdf=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
    if pc=='PM optical property':
        fdf=rdf[['compound',efcol]][rdf['pollutant_category']==pc][rdf[efcol].notna()].reset_index(drop=True)
        fdf=fdf.applymap(lambda x: rounding(x))
        return fdf
    else:
        fdf=rdf[['mm','formula','compound',efcol]][rdf['pollutant_category']==pc][rdf[efcol].notna()].reset_index(drop=True)
        fdf=fdf.applymap(lambda x: rounding(x))
        return fdf

# This function returns the EF of specified compound name and table name (integrated ef, processed ef and recommended ef)
def select_compound(ft, com_name,table_name):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')
    if table_name=='integrated ef':
        df=pd.read_sql(text('select * from Integrated_EF'), con=output_db)
        efcoldf=pd.read_sql(text('select * from bkdb_info_efcol'), con=bk_db)
        allcol= ['legend','fuel_type','measurement_type','MCE','EF']
        try:
                ind = get_ind (df, com_name)
                efcol=list(efcoldf['efcol'])
                efcoldf['EF']=df[efcol][df.index.isin(ind)].mean().values
                ll=efcoldf[allcol][efcoldf['fire_type']==ft]
                #ll=ll.sort_values(by='measurement_type')
                ll=ll[ll['EF'].notna()]
                ll=ll.sort_values(by='EF', ascending=False)
                ll=ll.reset_index(drop=True)
                ll=ll.applymap(lambda x: rounding(x))
                return ll
        except:
                return 'Compound not found. Search by formula'
    if table_name=='processed ef':
        df=pd.read_sql(text('select * from Processed_EF'), con=output_db)
        efcoldf=pd.read_sql(text('select * from info_efcol_processed_data'), con=bk_db)
        allcol= ['legend','fuel_type','measurement_type','MCE','EF']
        try:
            ind = get_ind (df, com_name)
            efcol=list(efcoldf['efcol'])
            efcoldf['EF']=df[efcol][df.index.isin(ind)].mean().values
            ll=efcoldf[allcol][efcoldf['fire_type']==ft]
            #ll=ll.sort_values(by='measurement_type')
            ll=ll[ll['EF'].notna()]
            ll=ll.sort_values(by='EF', ascending=False)
            ll=ll.reset_index(drop=True)
            ll=ll.applymap(lambda x: rounding(x))
            return ll
        except:
            return 'Compound not found. Search by formula'
    if table_name=='recommended ef':
        df=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
        df=df.applymap(lambda x: rounding(x))
        try:
            ind=get_ind_rdf (df,com_name)[0]
            col='AVG_'+ft.replace(' ','_')
            df[col]=df[col]
            return df[['mm','formula','compound',col]][ind:ind+1].reset_index(drop=True)
        except:
            return 'Compound not found. Search by formula'
    if table_name=='rdb':
        db_connection=connect_db('raw_db')
        tbl_info=pd.read_sql(text('select * from bkdb_info_rdb_ldb'), con=bk_db)
        tbl_info=tbl_info[tbl_info['fire_type']==ft][tbl_info['db']==table_name]
        tbl_ll=list(tbl_info['table'].unique())
        fdf=pd.DataFrame()
        try:
            for tbl in tbl_ll:
                rdf=pd.DataFrame()
                dd=pd.read_sql(text('select * from '+tbl), con=db_connection)
                efcol=list(tbl_info['efcol'][tbl_info['table']==tbl])
                ind=get_ind (dd, com_name)
                rdf['table']=[tbl]*len(efcol)
                rdf['efcol']=efcol
                rdf['EF']=dd[efcol].iloc[ind].mean().values
                rdf['db']=['raw db']*len(efcol)
                rdf=rdf.sort_values(by='EF', ascending=False).reset_index(drop=True)
                fdf=pd.concat([fdf,rdf])
            fdf=fdf[fdf['EF'].notna()].reset_index(drop=True)
            fdf=fdf.applymap(lambda x: rounding(x))
            return fdf
        except:
            return 'Compound not found. Search by formula'
    if table_name=='ldb':
        db_connection=connect_db('legacy_db')
        tbl_info=pd.read_sql(text('select * from bkdb_info_rdb_ldb'), con=bk_db)
        tbl_info=tbl_info[tbl_info['fire_type']==ft][tbl_info['db']==table_name]
        tbl_ll=list(tbl_info['table'].unique())
        fdf=pd.DataFrame()
        try:
            for tbl in tbl_ll:
                rdf=pd.DataFrame()
                dd=pd.read_sql(text('select * from '+tbl), con=db_connection)
                efcol=list(tbl_info['efcol'][tbl_info['table']==tbl])
                ind=get_ind (dd, com_name)
                rdf['table']=[tbl]*len(efcol)
                rdf['efcol']=efcol
                rdf['EF']=dd[efcol].iloc[ind].mean().values
                rdf['db']=['legacy db']*len(efcol)
                rdf=rdf.sort_values(by='EF', ascending=False).reset_index(drop=True)
                fdf=pd.concat([fdf,rdf])
            fdf=fdf[fdf['EF'].notna()].reset_index(drop=True)
            fdf=fdf.applymap(lambda x: rounding(x))
            return fdf
        except:
            return 'Compound not found. Search by formula'

        
def select_chemical_formula (ft, formula,table_name):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')
    if table_name=='integrated ef':
        df=pd.read_sql(text('select * from Integrated_EF'), con=output_db)
        efcoldf=pd.read_sql(text('select * from bkdb_info_efcol'), con=bk_db)
        ll=list(efcoldf['efcol'][efcoldf['fire_type']==ft])
        cols=['mm','formula','compound']+ll
        fdf=df[cols][df['formula']==formula].reset_index(drop=True)
        fdf=fdf.applymap(lambda x: rounding(x))
        return fdf
    if table_name=='processed ef':
        df=pd.read_sql(text('select * from Processed_EF'), con=output_db)
        efcoldf=pd.read_sql(text('select * from info_efcol_processed_data'), con=bk_db)
        ll=list(efcoldf['efcol'][efcoldf['fire_type']==ft])
        cols=['mm','formula','compound']+ll
        fdf=df[cols][df['formula']==formula].reset_index(drop=True)
        fdf=fdf.applymap(lambda x: rounding(x))
        return fdf
    if table_name=='recommended ef':
        output_db=connect_db('neiva_output_db')
        df=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
        df=df.applymap(lambda x: rounding(x))
        if isinstance(ft, list):
          cols = ['AVG_' + s.replace(' ','_') for s in ft]
        elif ft.lower()=='all':
          cols = df.columns[df.columns.str.contains('AVG_')]
        else:
          cols = ['AVG_'+ft.replace(' ','_')]        
        return df[['mm','formula', 'compound','id']+list(cols)][df['formula']==formula].reset_index(drop=True)
    

def abundant_nmog (ft, chem, aa):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')
    efcol='AVG_'+ft.replace(' ','_')
    rdf=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
    rdf=rdf[rdf[efcol].notna()]
    rdf=rdf[rdf['pollutant_category']=='NMOC_g'].reset_index(drop=True)
    pp=pd.read_sql(text('select * from Property_Surrogate'), con=output_db)
    rdf=rdf.sort_values(by=efcol, ascending=False)
    rdf=rdf.reset_index(drop=True)
    rdf=rdf.merge(pp[['id',chem, aa]], on='id', how='left')
    rdf=rdf.applymap(lambda x: rounding(x))
    return rdf[['mm','formula','compound', efcol, chem, aa, 'id']][:25].reset_index(drop=True)

def nmog_with_high_n (ft, chem, aa):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')
    ncol='N_'+ft.replace(' ','_')
    efcol='AVG_'+ft.replace(' ','_')
    rdf=pd.read_sql(text('select * from Recommended_EF'), con=output_db)
    rdf=rdf[rdf[efcol].notna()]
    rdf=rdf[rdf['pollutant_category']=='NMOC_g'].reset_index(drop=True)
    pp=pd.read_sql(text('select * from Property_Surrogate'), con=output_db)
    rdf=rdf.sort_values(by=ncol, ascending=False)
    rdf=rdf.reset_index(drop=True)
    rdf=rdf.merge(pp[['id',chem, aa]], on='id', how='left')
    rdf=rdf.applymap(lambda x: rounding(x))
    return rdf[['mm','formula','compound', efcol, ncol,chem, aa, 'id']][:25].reset_index(drop=True)

def ef_sorted_by_property (dd, ft, chem, model_surrogate, pr):
    output_db=connect_db('neiva_output_db')
    bk_db=connect_db('backend_db')
    nmog=join_ef_property(dd)
    # Set EF column based on the input parameter 'fire type'
    avgcol='AVG_'+ft.replace(' ','_')
    ncol='N_'+ft.replace(' ','_')
    stdcol='STD_'+ft.replace(' ','_')
    nmog['ef']=nmog[avgcol]
    nmog=nmog[nmog['ef'].notna()].reset_index(drop=True)
    nmog['mole']=nmog['ef']/nmog['mm']
    nmog['mole_frac']=nmog['mole']/nmog['mole'].sum()
    nmog=nmog[nmog[chem]==model_surrogate]
    nmog=nmog.sort_values(by=pr, ascending=False)
    nmog=nmog.reset_index(drop=True)
    nmog=nmog.applymap(lambda x: rounding(x))
    nmog=nmog[['mm','formula','compound',avgcol,ncol,stdcol, chem, pr]]
    return nmog
    
def compare_lab_field (ft, com_name,table_name):
    bk_db=connect_db('backend_db')
    output_db=connect_db('neiva_output_db')

    if table_name=='integrated ef':
        df=pd.read_sql(text('select * from Integrated_EF'), con=output_db)
        efcoldf=pd.read_sql(text('select * from bkdb_info_efcol'), con=bk_db)

    if table_name=='processed ef':
        df=pd.read_sql(text('select * from Processed_EF'), con=output_db)
        efcoldf=pd.read_sql(text('select * from info_efcol_processed_data'), con=bk_db)
    try:
        iind = get_ind (df, com_name)
        
        co_ind=df[df['id']=='InChI=1S/CO/c1-2'].index[0]
        co2_ind=df[df['id']=='InChI=1S/CO2/c2-1-3'].index[0]
        
        efcol=list(efcoldf['efcol'])

        efcoldf[com_name]=df[efcol][df.index.isin(iind)].mean().values
        efcoldf['co']=df[efcol].iloc[co_ind].values
        efcoldf['co2']=df[efcol].iloc[co2_ind].values
        
        efcoldf['MCE']=(efcoldf['co2']/44) / ((efcoldf['co2']/44) + (efcoldf['co']/28))
    
        ll=efcoldf[[com_name,'MCE','measurement_type']][efcoldf['fire_type']==ft]
        ll=ll.sort_values(by='measurement_type')
        ll=ll[ll[com_name].notna()]
        ll=ll.reset_index(drop=True)
        
        lab_avg=ll[com_name][ll['measurement_type']=='lab'].mean()
        field_avg=ll[com_name][ll['measurement_type']=='field'].mean()
        mce_lab = ll['MCE'][ll['measurement_type']=='lab'].mean()
        mce_field = ll['MCE'][ll['measurement_type']=='field'].mean()
        n_lab=len(ll[com_name][ll['measurement_type']=='lab'])
        n_field=len(ll[com_name][ll['measurement_type']=='field'])
        
        data = {
            'Mean': [com_name.capitalize()+' EF', 'MCE', 'data count'],
            'Lab': [lab_avg, mce_lab, n_lab],
            'Field': [field_avg,mce_field,n_field]
             }          
        df=pd.DataFrame(data)
        df=df.applymap(lambda x: rounding(x))
        return df
    except:
        return 'Compound not found. Search by formula'

    


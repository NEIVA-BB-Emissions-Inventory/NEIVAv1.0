#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 14:51:59 2024

@author: samiha
"""

import pandas as pd
import numpy as np
import pubchempy as pcp


def assign_id (df):
    for i in range(len(df)):
        try:
            c=pcp.get_compounds(df['compound'].iloc[i], 'name')
            df.loc[i,'id']=c[0].inchi
        except:
            df.loc[i,'id']=np.nan
    return df

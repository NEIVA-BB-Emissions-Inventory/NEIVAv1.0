#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 27 20:34:14 2022

@author: samiha
"""

from prettytable import PrettyTable
import pandas as pd
import numpy as np


def create_PrettyTable_col2(list_col,col_list1,col_list2):
    t = PrettyTable(list_col)
    allrow=[]
    for k in range(len(col_list1)):
        single_rowlist=[]
        single_rowlist.append(col_list1[k])
        single_rowlist.append(col_list2[k])
        allrow.append(single_rowlist)
    for i in range(len(allrow)):
        t.add_row(allrow[i])
    return t
    
def create_PrettyTable_col3(list_col,col_list1,col_list2,col_list3):
    t = PrettyTable(list_col)
    allrow=[]
    for k in range(len(col_list1)):
        single_rowlist=[]
        single_rowlist.append(col_list1[k])
        single_rowlist.append(col_list2[k])
        single_rowlist.append(col_list3[k])
        allrow.append(single_rowlist)
    for i in range(len(allrow)):
        t.add_row(allrow[i])
    return t

def create_PrettyTable_col4(list_col,col_list1,col_list2,col_list3,col_list4):
    t = PrettyTable(list_col)
    allrow=[]
    for k in range(len(col_list1)):
        single_rowlist=[]
        single_rowlist.append(col_list1[k])
        single_rowlist.append(col_list2[k])
        single_rowlist.append(col_list3[k])
        single_rowlist.append(col_list4[k])
        allrow.append(single_rowlist)
    for i in range(len(allrow)):
        t.add_row(allrow[i])
    return t



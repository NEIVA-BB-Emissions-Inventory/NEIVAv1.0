#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 14:20:03 2023

@author: samiha
"""

# import pandas as pd
from sqlalchemy import create_engine
import mysql.connector
import pandas as pd

def connect_db(db):
    """
    Connects to a MySQL database using SQLAlchemy.
    Input:
    - db (str): Database name (e.g., legacy_db, raw_db, primary_db).
    
    Returns:
        The connection. 
    """
    engine = create_engine("mysql+pymysql://root:root@localhost/"+db)
    con = engine.connect()
    return con

def get_table_name(database_name):
    '''
    Connects to the specified MySQL database and creates a cursor for queries.
    
    Input:
        database_name (str): Name of the target database.
    
    Returns:
        Cursor object for executing queries in the database.
    '''    
    # Define your connection parameters
    # config = {
    #     'user': 'root',
    #     'password': 'Samiha_056',
    #     'host': 'localhost',  # adjust if your MySQL server isn't on localhost
    #     'database': database_name,
    # }
    # Establish the connection
    conn = mysql.connector.connect(user='root', password='root', host='localhost', database=database_name)
    # conn = mysql.connector.connect(**config)
    
    # Create a cursor object using the connection
    cursor = conn.cursor()
    
    query = "SHOW TABLES"
    cursor.execute(query)
    
    # Fetch and print table names
    tables = cursor.fetchall()
    
    table_name=[]
    for tbl in tables:
        # table_name = pd.concat([table_name,tbl[0]], ignore_index=True)
        table_name.append(tbl[0])
    
    # Close the cursor and connection
    cursor.close()
    conn.close()
    # Create and return a cursor object using the connection
    return table_name












# -*- coding: utf-8 -*-
"""
Created on Thu May 14 11:18:21 2020

@author: Swapnil Keshari

Summary: This file reads the data form the csv file and operates on it to create
another dataframe which has unique rows and Phospho as columns
"""

import pandas as pd

ph = pd.read_csv('D:\Desktop\IITB\Project_Juilee\phosphoproteomics.csv', nrows=50)

#Replacing all the NA with '0'
ph.fillna(0, inplace=True)  
ph.sort_values(by=['Gene','Sequence'],inplace=True)
ph.reset_index(drop=True, inplace=True)
ph.drop(columns = "Unnamed: 0", inplace = True)

#Modifying the 'Modification' and extracting Phophorylation sites from it using string manipulations
new= ph["Modifications"].str.split(' (?=\dxPhospho)', n = 1, expand = True)
ph["Phospho"] = new[1]
ph['Phospho'] = ph["Phospho"].astype(str)
ph['Phospho'] = ph["Phospho"].str.replace('\[|\([^)]*\)|\]|Phospho\s', '').str.strip()
ph['Phospho'] = ph["Phospho"].str.replace('\;\s', '-')
ph['Phospho'] = ph["Phospho"].str.replace('T/S', 'S/T')
ph.drop(columns =["Modifications"], inplace = True)

#Creating a Sequence1 which is a unique identifier for each modification
ph['Sequence1']=ph['Gene'].astype(str)+"_"+ph['Sequence'].astype(str)+"_"+ph['Phospho'].astype(str)
#Grouping all the non unique entries and adding the values corresponding to it    
res = ph.groupby('Sequence1', as_index=False).sum()
#Creating  new dataframe having separate columns for Gene Sequence and Phospho
new= res["Sequence1"].str.split('_', n = 2, expand = True)
res["Gene"] = new[0].astype(str)
res["Sequence"]=new[1].astype(str)
res["Phospho"]=new[2].astype(str)

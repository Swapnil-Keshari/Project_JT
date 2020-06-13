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

#ph.dropna(inplace = True) 
#Modifying the 'Modification' and extracting Phophorylation sites from it using string manipulations
new= ph["Modifications"].str.split(' (?=\dxPhospho)', n = 1, expand = True)
ph["Phospho"] = new[1]
ph['Phospho'] = ph["Phospho"].astype(str)
ph['Phospho'] = ph["Phospho"].str.replace('\[|\([^)]*\)|\]|Phospho\s', '').str.strip()
ph['Phospho'] = ph["Phospho"].str.replace('\;\s', '-')
ph['Phospho'] = ph["Phospho"].str.replace('T/S', 'S/T')
ph.drop(columns =["Modifications"], inplace = True)

#ph['Sequence1'] = ph[['Gene', 'Sequence', 'Phospho']].agg('-'.join, axis=1)
#Creating a Sequence1 which is a unique identifier for each modification
ph['Sequence1']=ph['Gene'].astype(str)+"_"+ph['Sequence'].astype(str)+"_"+ph['Phospho'].astype(str)
#Grouping all the non unique entries and adding the values corresponding to it    
res = ph.groupby('Sequence1', as_index=False).sum()
#Creating  new dataframe having separate columns for Gene Sequence and Phospho
new= res["Sequence1"].str.split('_', n = 2, expand = True)
res["Gene"] = new[0].astype(str)
res["Sequence"]=new[1].astype(str)
res["Phospho"]=new[2].astype(str)
























#new= ph["Modifications"].str.split("dation", n = 1, expand = True)
#new= new[1].str.split("];", n = 1, expand = True)
#ph["Oxidation"] = new[0]
#ph['Oxidation'] = ph["Oxidation"].str.replace('\[|\]\;', '').str.strip()
#ph['Oxidation'] = ph["Oxidation"].str.replace('\;\s', '-')

#new= ph["Modifications"].str.split("methyl", n = 1, expand = True) 
#new= new[1].str.split("];", n = 1, expand = True)
#ph["Carbamidomethyl"] = new[0]
#ph["Carbamidomethyl"] = ph["Carbamidomethyl"].str.replace('\[|\]\;', '').str.strip()
#ph["Carbamidomethyl"] = ph["Carbamidomethyl"].str.replace('\;\s', '-')

#new= ph["Modifications"].str.split("pho", n = 1, expand = True) 
#ph["Phospho"] = new[1]
#ph['Phospho'] = ph["Phospho"].astype(str)
#ph['Phospho'] = ph["Phospho"].str.replace('\[|\([^)]*\)|\]', '').str.strip()
#ph['Phospho'] = ph["Phospho"].str.replace('\;\s', '-')
#ph['Phospho'] = ph["Phospho"].str.replace('T/S', 'S/T')

#ph.drop(columns =["Modifications"], inplace = True)

#for index in range(len(ph['Sequence'])):
#    ph['Sequence'][index]=str(ph['Gene'][index])+"_"+str(ph['Sequence'][index])+"_"+str(ph['Phospho'][index])
#    +"_"+str(ph['Oxidation'][index])+"_" +str(ph['Carbamidomethyl'][index])

#ph['Phospho'] = ph["Phospho"].astype(str)
##ph['Phospho'] = ph["Phospho"].str.split('[', 1).str[1].str.strip()
##ph['Phospho'] = ph["Phospho"].str.split('(', 1).str[0].str.strip()
##ph['Phospho'] = ph["Phospho"].str.split(']', 1).str[0].str.strip()
#ph['Phospho'] = ph["Phospho"].str.replace('\[|\([^)]*\)|\]', '').str.strip()
#ph['Phospho'] = ph["Phospho"].str.replace('\;\s', '-')

##Creation of dictionary
#s= pd.Series(ph["Gene"]) 
#Dict={}
#Dict={key:[] for key in s.unique()}
#print len(Dict.keys())
#iteration=1
#for key in Dict.keys():
#    print iteration
#    iteration+=1
#    for index in range(len(ph["Gene"])):
#        if key == ph["Gene"][index]:
#            Dict[key].append(str(ph['Oxidation'][index])+"_" +str(ph['Carbamidomethyl'][index])+"_"+str(ph['Phospho'][index]))
            

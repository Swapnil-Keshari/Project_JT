# -*- coding: utf-8 -*-
"""
Created on Thu May 14 11:18:21 2020

@author: 98kes
"""

import pandas as pd

ph = pd.read_excel('D:\Desktop\IITB\Project_Juilee\phosphoproteomics.xlsx',nrows=12)

ph.fillna(0, inplace=True)  
ph.sort_values(by=['Gene','Sequence'],inplace=True)
ph.reset_index(drop=True, inplace=True)
ph.drop(columns = "Unnamed: 0", inplace = True)
#ph.dropna(inplace = True) 
#ph.reset_index(drop=True, inplace=True)

new= ph["Modifications"].str.split("dation", n = 1, expand = True)
new= new[1].str.split("];", n = 1, expand = True)
ph["Oxidation"] = new[0]
ph['Oxidation'] = ph["Oxidation"].str.replace('\[|\]\;', '').str.strip()
ph['Oxidation'] = ph["Oxidation"].str.replace('\;\s', '-')

#new= ph["Modifications"].str.split("methyl", n = 1, expand = True) 
#new= new[1].str.split("];", n = 1, expand = True)
#ph["Carbamidomethyl"] = new[0]
#ph["Carbamidomethyl"] = ph["Carbamidomethyl"].str.replace('\[|\]\;', '').str.strip()
#ph["Carbamidomethyl"] = ph["Carbamidomethyl"].str.replace('\;\s', '-')

new= ph["Modifications"].str.split("pho", n = 1, expand = True) 
ph["Phospho"] = new[1]
ph['Phospho'] = ph["Phospho"].astype(str)
ph['Phospho'] = ph["Phospho"].str.replace('\[|\([^)]*\)|\]', '').str.strip()
ph['Phospho'] = ph["Phospho"].str.replace('\;\s', '-')
ph['Phospho'] = ph["Phospho"].str.replace('T/S', 'S/T')

ph.drop(columns =["Modifications"], inplace = True)

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
            

# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 15:13:10 2020

@author: Swapnil Keshari
"""
import numpy as np
import pandas as pd
from dataprocessing import ph

ph['Sequence1']=ph['Sequence']

for index in range(len(ph['Sequence'])):
    ph['Sequence1'][index]=str(ph['Gene'][index])+"_"+str(ph['Sequence'][index])+"_"+str(ph['Phospho'][index])+"_"+str(ph['Oxidation'][index])

#Creation of dictionary
s= pd.Series(ph["Gene"]) 
Dict={}
Dict_1={}
Dict={gene:[] for gene in s.unique()}
print len(Dict.keys())
for gene in Dict.keys():
    for index in range(len(ph["Gene"])):
        if gene == ph["Gene"][index]:
            Dict[gene].append(ph['Sequence'][index])
    Dict[gene]=list(set(Dict[gene]))
    Dict_1={seq:[] for seq in Dict[gene]}
    for seq in Dict_1.keys():
        for index in range(len(ph["Gene"])):
            if seq == ph["Sequence"][index]:
                Dict_1[seq].append(ph['Phospho'][index])
        Dict_1[seq]=list(set(Dict_1[seq]))
        for st in Dict_1[seq]:
            flag=False
            if 'S/T' in st:
                if '-' in st:
                    s=st.replace('\-S/T-','--')
                    s=st.replace('\-S/T|\S/T-','-')
                    if s in Dict_1[seq]:
                        flag=True
            if flag:
                Dict_1[seq].remove(st)
        
                        

        
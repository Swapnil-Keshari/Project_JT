# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 15:13:10 2020

@author: Swapnil Keshari
"""
import numpy as np
import pandas as pd
import re
from dataprocessing import ph

ph['Sequence1']=ph['Sequence']

for index in range(len(ph['Sequence'])):
    ph['Sequence1'][index]=str(ph['Gene'][index])+"_"+str(ph['Sequence'][index])+"_"+str(ph['Phospho'][index])

res = ph.groupby('Sequence1', as_index=False).sum()
#Creation of dictionary
Dict={}
Dict={gene:{} for gene in set(ph['Gene'])}
print len(Dict.keys())
for gene in Dict.keys():
    for index in range(len(ph["Gene"])):
        if gene == ph["Gene"][index]:
            Dict[gene][ph['Sequence'][index]]=[]
    for seq in Dict[gene].keys():
        for index in range(len(ph["Gene"])):
            if seq == ph["Sequence"][index]:
                Dict[gene][seq].append(ph['Phospho'][index])
        for ele in range(len(Dict[gene][seq])):
            if (re.match('1xS\d|1xT\d|1xS/T',Dict[gene][seq][ele]) is None) and '1x' in Dict[gene][seq][ele]:
                for ele1 in range(len(Dict[gene][seq])):
                    if (re.match('1xS\d|1xT\d',Dict[gene][seq][ele1]) is not None) and '1x' in Dict[gene][seq][ele1]:
                        Dict[gene][seq][ele] = Dict[gene][seq][ele1]
            elif 'S/T' in Dict[gene][seq][ele] and '2x' in Dict[gene][seq][ele] and '2xS/T' not in Dict[gene][seq][ele]:
                s = Dict[gene][seq][ele].replace('2x','')
                for value in s.split('-'):
                    if (re.match('S\d|T\d',value) is not None):
                        ele2 =value
                    for ele1 in range(len(Dict[gene][seq])):
                        if ele2 in Dict[gene][seq][ele1] and '2x' in Dict[gene][seq][ele1] and 'S/T' not in Dict[gene][seq][ele1]:
                            Dict[gene][seq][ele] = Dict[gene][seq][ele1]

#                m = re.search('S\d\d', ele)
#                if m:
#                    found = m.group(1)
#                    print found
                    


#    Dict[gene]=list(set(Dict[gene]))

#Dict_1={seq:[] for seq in Dict[gene]}
#    for seq in Dict_1.keys():
#        for index in range(len(ph["Gene"])):
#            if seq == ph["Sequence"][index]:
#                Dict_1[seq].append(ph['Phospho'][index])
#        Dict_1[seq]=list(set(Dict_1[seq]))
#        for st in Dict_1[seq]:
#            flag=False
#            if 'S/T' in st:
#                if '-' in st:
#                    s=st.replace('\-S/T-','--')
#                    s=st.replace('\-S/T|\S/T-','-')
#                    if s in Dict_1[seq]:
#                        flag=True
#            if flag:
#                Dict_1[seq].remove(st)   











 
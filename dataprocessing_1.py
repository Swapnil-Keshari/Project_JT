# -*- coding: utf-8 -*-
"""
Created on Fri May 29 08:28:51 2020

@author: Swapnil Keshari

Summary: This file imports the res dataframe and operates on it to delete the unknown phosphorylation
sites and add it to the known site by matching. Finally the result is a multi layer dictionary
Dict2={Gene:{Sequence:{Phospho:[[dataframe]]}}}
"""

import pandas as pd
import numpy as np
import re
from dataprocessing import res

#Creating dictionary from the initial dataframe imported
Dict1=res.to_dict('series')
#This is the start of multilevel dictionary in which genes are mapped to their sequences which are mapped to the corresponding
#phosphorylation sites. The sites have the initial values stored corresponding to them in form of data frame
Dict2={}
Dict2={gene:{} for gene in set(res['Gene'])}
#Adds sequences corresponding to the gene
for index in range(len(Dict1["Gene"])):
    print index
    if Dict1["Gene"][index] in Dict2.keys():
        Dict2[Dict1["Gene"][index]][Dict1['Sequence'][index]]={}
#Adds Phospo sites corresponding to the sequence and then appending the initial values as dataframe
for index in range(len(Dict1["Gene"])):
    print str (index) + "A"
    if Dict1["Gene"][index] in Dict2.keys():
        if Dict1["Sequence"][index] in Dict2[Dict1["Gene"][index]].keys():
            Dict2[Dict1["Gene"][index]][Dict1['Sequence'][index]][Dict1['Phospho'][index]]=[]
            Dict2[Dict1["Gene"][index]][Dict1['Sequence'][index]][Dict1['Phospho'][index]].append(res.iloc[[index]].select_dtypes(exclude=[np.object]))
#This loop looks for all the unknown phosphorylation sites and add it to similar one and deletes them
for gene in Dict2.keys():
    print gene
    for seq in Dict2[gene].keys():
        print seq
        for ele in Dict2[gene][seq].keys():
            #For 1xS or 1xT type of entries. We don't touch 1xS/T
            if (re.match('1xS\d|1xT\d|1xS/T', ele) is None) and '1x' in ele: #Check for 1xS,1xT AND 1x
                for ele1 in Dict2[gene][seq].keys():
                    if (re.match('1xS\d|1xT\d',str(ele1)) is not None) and '1x' in str(ele1):# Check for 1xS<digit>,1xT<digit> AND 1x
                        #We add the values of the dataframe to the first match and then delete the element ele
                        a= pd.concat([Dict2[gene][seq][ele1][0],Dict2[gene][seq][ele][0]]).groupby(level=0).sum()
                        Dict2[gene][seq][ele1][0]= a.reset_index().replace({'index':{Dict2[gene][seq][ele][0].index.values[0]:Dict2[gene][seq][ele1][0].index.values[0]}}).groupby('index').sum()
                        del Dict2[gene][seq][ele]
                        break
            #For 2xS/T-S<digit> or 2xS/T-T<digit> type of entries. We don't touch 2xS/T
            if 'S/T' in ele and '2x' in ele and '2xS/T' not in ele:# Check for 2xS/T-S<digit> or 2xS/T-T<digit> type
                s = ele.replace('2x','')
                for value in s.split('-'):
                    if (re.match('S\d|T\d',value) is not None):# Extracting S<digit> or T<digit>
                        ele2 =value
                for ele1 in Dict2[gene][seq].keys():
                    if ele2 in ele1 and '2x' in ele1 and 'S/T' not in ele1:
                        #We add the values of the dataframe to the first match and then delete the element ele
                        a= pd.concat([Dict2[gene][seq][ele1][0],Dict2[gene][seq][ele][0]]).groupby(level=0).sum()
                        Dict2[gene][seq][ele1][0]= a.reset_index().replace({'index':{Dict2[gene][seq][ele][0].index.values[0]:Dict2[gene][seq][ele1][0].index.values[0]}}).groupby('index').sum()
                        del Dict2[gene][seq][ele]
                        break
            #For 3x which has S/T in it type of entries. We don't touch 3xS/T
            elif 'S/T' in ele and '3x' in ele and '3xS/T' not in ele:
                t = ele.replace('3x','')
#                print t
                u=t.split('-')
                u.remove('S/T')#remove S/T from the set so that only digit values remain
                v=set(t)
#                print v
                for ele1 in Dict2[gene][seq].keys():
                    if '3x' in ele1 and 'S/T' not in ele1:
                        t1 = ele1.replace('3x','')
                        v1 = set (t1.split('-'))
                        if v.issubset(v1):# if the digit values are available in the other 3x
                            #We add the values of the dataframe to the first match and then delete the element ele
                            a= pd.concat([Dict2[gene][seq][ele1][0],Dict2[gene][seq][ele][0]]).groupby(level=0).sum()
                            Dict2[gene][seq][ele1][0]= a.reset_index().replace({'index':{Dict2[gene][seq][ele][0].index.values[0]:Dict2[gene][seq][ele1][0].index.values[0]}}).groupby('index').sum()
                            del Dict2[gene][seq][ele]
                            break

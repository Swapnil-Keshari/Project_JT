# -*- coding: utf-8 -*-
"""
Created on Fri May 29 08:28:51 2020

@author: sanka
"""

import pandas as pd
import numpy as np
import re
from dataprocessing import res

#ph['Sequence1']=ph['Sequence']
#
#for index in range(len(ph['Sequence'])):
#    ph['Sequence1'][index]=str(ph['Sequence'][index])+"_"+str(ph['Phospho'][index])
#
#res = ph.groupby('Sequence1', as_index=False).sum()
Dict1=res.to_dict('series')

Dict2={}
Dict2={gene:{} for gene in set(res['Gene'])}

#Dict3={}
#Dict3={gene:{} for gene in set(res['Gene'])}
for index in range(len(Dict1["Gene"])):
    print index
    if Dict1["Gene"][index] in Dict2.keys():
        Dict2[Dict1["Gene"][index]][Dict1['Sequence'][index]]={}

for index in range(len(Dict1["Gene"])):
    print str (index) + "A"
    if Dict1["Gene"][index] in Dict2.keys():
        if Dict1["Sequence"][index] in Dict2[Dict1["Gene"][index]].keys():
            Dict2[Dict1["Gene"][index]][Dict1['Sequence'][index]][Dict1['Phospho'][index]]=[]
            Dict2[Dict1["Gene"][index]][Dict1['Sequence'][index]][Dict1['Phospho'][index]].append(res.iloc[[index]].select_dtypes(exclude=[np.object]))

for gene in Dict2.keys():
    print gene
#    for index in range(len(Dict1["Gene"])):
#        print index
#        if gene == Dict1["Gene"][index]:
#            Dict2[gene][Dict1['Sequence'][index]]={}
#    for seq in Dict2[gene].keys():
#        print 'Seq'
#        for index in range(len(Dict1["Gene"])):
#            if seq == Dict1["Sequence"][index]:
#                Dict2[gene][seq][Dict1['Phospho'][index]]=[]
#        for phospho in Dict2[gene][seq]:
#            print 'Phospho'
#            for index in range(len(Dict1["Gene"])):
#                if (str(gene)+'_'+str(seq)+'_'+str(phospho)) == Dict1["Sequence1"][index]:
#                    Dict2[gene][seq][phospho].append(res.iloc[[index]].select_dtypes(exclude=[np.object]))
#    print 'Dict Formed'
    for seq in Dict2[gene].keys():
        print seq
        for ele in Dict2[gene][seq].keys():
#            if len(Dict2[gene][seq][ele])>1:
#                x = len(Dict2[gene][seq][ele])-1
#                for i in range(x):
#                    Dict2[gene][seq][ele][0]= pd.concat([Dict2[gene][seq][ele][0],Dict2[gene][seq][ele][1]]).groupby(level=0).sum()
#                    Dict2[gene][seq][ele][0]= Dict2[gene][seq][ele][0].reset_index().replace({'index':{Dict2[gene][seq][ele][0].index.values[0]:Dict2[gene][seq][ele][1].index.values[0]}}).groupby('index').sum()
#                    del Dict2[gene][seq][ele][1]
                
            if (re.match('1xS\d|1xT\d|1xS/T', ele) is None) and '1x' in ele: #S/T
                for ele1 in Dict2[gene][seq].keys():
                    if (re.match('1xS\d|1xT\d',str(ele1)) is not None) and '1x' in str(ele1):
                        a= pd.concat([Dict2[gene][seq][ele1][0],Dict2[gene][seq][ele][0]]).groupby(level=0).sum()
                        Dict2[gene][seq][ele1][0]= a.reset_index().replace({'index':{Dict2[gene][seq][ele][0].index.values[0]:Dict2[gene][seq][ele1][0].index.values[0]}}).groupby('index').sum()
                        del Dict2[gene][seq][ele]
                        break
            if 'S/T' in ele and '2x' in ele and '2xS/T' not in ele:
                s = ele.replace('2x','')
                for value in s.split('-'):
                    if (re.match('S\d|T\d',value) is not None):
                        ele2 =value
                for ele1 in Dict2[gene][seq].keys():
                    if ele2 in ele1 and '2x' in ele1 and 'S/T' not in ele1:
                        a= pd.concat([Dict2[gene][seq][ele1][0],Dict2[gene][seq][ele][0]]).groupby(level=0).sum()
                        Dict2[gene][seq][ele1][0]= a.reset_index().replace({'index':{Dict2[gene][seq][ele][0].index.values[0]:Dict2[gene][seq][ele1][0].index.values[0]}}).groupby('index').sum()
                        del Dict2[gene][seq][ele]
                        break
            elif 'S/T' in ele and '3x' in ele and '3xS/T' not in ele:
                t = ele.replace('3x','')
#                print t
                u=t.split('-')
                u.remove('S/T')
                v=set(t)
#                print v
                for ele1 in Dict2[gene][seq].keys():
                    if '3x' in ele1 and 'S/T' not in ele1:
                        t1 = ele1.replace('3x','')
                        v1 = set (t1.split('-'))
                        if v.issubset(v1):
                            a= pd.concat([Dict2[gene][seq][ele1][0],Dict2[gene][seq][ele][0]]).groupby(level=0).sum()
                            Dict2[gene][seq][ele1][0]= a.reset_index().replace({'index':{Dict2[gene][seq][ele][0].index.values[0]:Dict2[gene][seq][ele1][0].index.values[0]}}).groupby('index').sum()
                            del Dict2[gene][seq][ele]
                            break
                

#df_add = df1.add(df2, fill_value=0)


























#from numpy.random import randint
#
#df = pd.DataFrame(columns=['lib','lib1', 'qty1', 'qty2','qty3','qty4','qty5'])
#for i in range(5):
#     df.loc[i] = ['name']+ ['nam'+str(i)] + list(randint(10, size=5))
#ph['Zeroes']=(ph == 0).astype(int).sum(axis=1)
#ph_1 = ph[ph['Zeroes'] > 20] 
#res = ph_1.groupby('Sequence', as_index=False).sum()
#res.drop(columns =["Zeroes"], inplace = True)
#res['Zeroes']=(res == 0).astype(int).sum(axis=1)
##rslt_df = res.loc[res['Zeroes'] > 20] 
#res1 =ph_1.groupby('Sequence').agg('sum')
#print ph
#print res1
#print res
#df = pd.DataFrame({'PC': ['DE101','DE101'], 'RatingCY': [None,"AA+"], 'RatingPY': ['AA', None],'HT':['GV','GV'],'MV1':[5.0,3.0],'MV2':[1.0,2.0]}, )
#
#a=df.head(1).combine_first(df.tail(1))
#
#obj_df = df.select_dtypes(include=[np.object])
#num_df = df.select_dtypes(exclude=[np.object])
#
#b= obj_df.head(1).combine_first(obj_df.tail(1)).join(num_df.head(1)+(num_df.tail(1)))
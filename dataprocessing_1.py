# -*- coding: utf-8 -*-
"""
Created on Fri May 29 08:28:51 2020

@author: sanka
"""

import pandas as pd
from dataprocessing import ph
#from numpy.random import randint
#
#df = pd.DataFrame(columns=['lib','lib1', 'qty1', 'qty2','qty3','qty4','qty5'])
#for i in range(5):
#     df.loc[i] = ['name']+ ['nam'+str(i)] + list(randint(10, size=5))
ph['Zeroes']=(ph == 0).astype(int).sum(axis=1)
ph_1 = ph[ph['Zeroes'] > 20] 
res = ph_1.groupby('Sequence', as_index=False).sum()
res.drop(columns =["Zeroes"], inplace = True)
res['Zeroes']=(res == 0).astype(int).sum(axis=1)
#rslt_df = res.loc[res['Zeroes'] > 20] 
res1 =ph_1.groupby('Sequence').agg('sum')
print ph
print res1
print res
#df = pd.DataFrame({'PC': ['DE101','DE101'], 'RatingCY': [None,"AA+"], 'RatingPY': ['AA', None],'HT':['GV','GV'],'MV1':[5.0,3.0],'MV2':[1.0,2.0]}, )
#
#a=df.head(1).combine_first(df.tail(1))
#
#obj_df = df.select_dtypes(include=[np.object])
#num_df = df.select_dtypes(exclude=[np.object])
#
#b= obj_df.head(1).combine_first(obj_df.tail(1)).join(num_df.head(1)+(num_df.tail(1)))
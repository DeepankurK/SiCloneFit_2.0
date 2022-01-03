#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 14:24:06 2020

@author: deepank
"""

import scipy.special as sp
import numpy as np  
import math
mut='/home/deepank/Downloads/Prof_Hamim/Data/P342/P342_Mut.txt'
ref='/home/deepank/Downloads/Prof_Hamim/Data/P342/P342_Ref.txt'

dst=open(mut,"rt")
m=dst.readlines()
dst.close()
dst=open(ref,"rt")
r=dst.readlines()
dst.close()
c=[]
q=np.zeros((1430,16),dtype=np.float)
for i in range(1,1431):
    c.append(r[i].split('\n')[0].split('\t')[0])
    for j in range(1,17):
        a=int(m[i].split('\n')[0].split('\t')[j])
        b=int(r[i].split('\n')[0].split('\t')[j])
        d=0
        for k in range(min(a,b)+1,a+b+1):
            d=d+math.log(k)-math.log(a+b+1-k)
        if (d<0.0001): d=0
        q[i-1][j-1]=d
import pandas as pd 
df=pd.DataFrame(q)
df['cell']=c
df.to_csv('P342_log_try.csv')

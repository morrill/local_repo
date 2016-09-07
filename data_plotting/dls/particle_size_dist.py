# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 15:21:27 2016

@author: Michael Orrill
"""
#%%
import numpy as np
import pandas as pd
from scipy import stats, integrate
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(color_codes=True)
#%%
st1 = np.loadtxt('standard1.asc', delimiter="\t", skiprows=89)
st2 = np.loadtxt('standard2.asc', delimiter="\t", skiprows=89)

#%% Trying to plot it with a bar plot

plt.figure(figsize=(45,10))
f, (ax1, ax2) = plt.subplots(2)
sns.barplot(st1[:,0],st1[:,1], ax=ax1, palette="PuBuGn_d")
plt.xticks(rotation=90)
plt.xlabel('Particle Size [nm]')
plt.ylabel('Relative Intensity')
sns.barplot(st2[:,0],st2[:,1], ax=ax2, palette="YlGnBu")
plt.xticks(rotation=90)
plt.xlabel('Particle Size [nm]')
plt.ylabel('Relative Intensity')
plt.show()
#%% Try and change the data so that you can make a histogram. Make 10,000 the highest number of particles counted per one size. 
# then take intensity*10,000 and add that many of the size number to an array.

intensity = Sonicatedcsv[:,1]
size = Sonicatedcsv[:,0]
num = 100
datapool =[]

for i in range(len(size)):
    temp_num = intensity[i]*num
    if temp_num != 0:
        for k in range(int(temp_num)):
            datapool.append(temp_num)
            
sns.distplot(datapool)
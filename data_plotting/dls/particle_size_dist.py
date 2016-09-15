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

#%% read in polystyrene standards from 9/6/16
st3 = np.loadtxt('standard3_9_6_16.asc', delimiter="\t", skiprows=89)
st4 = np.loadtxt('standard4_9_6_16.asc', delimiter="\t", skiprows=89)
st5 = np.loadtxt('standard5_9_6_16_new_water.asc', delimiter="\t", skiprows=89)
st6 = np.loadtxt('standard6_9_6_16_new_water.asc', delimiter="\t", skiprows=89)

#%% read in carbon ink from 9/1/16
c41 = np.loadtxt('C4.1_30ul_9_1_16.asc', delimiter="\t", skiprows=89)
c42 = np.loadtxt('C4.1_30ul_9_1_16.asc', delimiter="\t", skiprows=89)
c43 = np.loadtxt('C4.1_30ul_9_1_16.asc', delimiter="\t", skiprows=89)

#%% read in carbon ink from 9/6/16
c44 = np.loadtxt('C4.1_sonicated_9_6_16.asc', delimiter="\t", skiprows=89)

#%% Plot standards with a bar plot
plt.figure(figsize=(15,5))
plt.bar(st3[:,0],st3[:,1], width=15, color=(0,0.6,0.2,0.2))
plt.bar(st4[:,0],st4[:,1], width=15, color=(0,1,1,0.2))
plt.bar(st5[:,0],st5[:,1], width=15, color=(0,0.1,0.5,0.2))
plt.bar(st6[:,0],st6[:,1], width=15, color=(0.8,0.5,1,0.2))
plt.xlabel('Particle Size [nm]')
plt.ylabel('Relative Intensity')
plt.xlim([0,800])
plt.title('Polystyrene standard Size = 495 nm')
plt.grid()

#%% plot carbon data from 9/1/16


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
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 19:17:58 2016

@author: Michael Orrill

This is for plotting results from rheometer viscosity measurements
"""
import numpy as np
import matplotlib.pyplot as plt
#%%
rt5 = np.loadtxt('RT5_standard.txt', skiprows=3, usecols=(0,1,2))
#eg11 = np.loadtxt('PureEG1_1', skiprows=3, usecols=(0,1,2)) # pure ethylene glycol 1.1
#eg12 = np.loadtxt('PureEG1_2', skiprows=3, usecols=(0,1,2)) # pure ethylene glycol 1.1
#eg21 = np.loadtxt('PureEG2_1', skiprows=3, usecols=(0,1,2)) # pure ethylene glycol 1.1
#eg22 = np.loadtxt('PureEG2_2', skiprows=3, usecols=(0,1,2)) # pure ethylene glycol 1.1

#%% plot carbon data from 9/1/16
plt.figure(figsize=(15,5))
plt.bar(st_auto_on2_200[:,0],st_auto_on2_200[:,1], width=15, color=(0.8,0.5,1,0.2))
plt.xlabel('Particle Size [nm]')
plt.ylabel('Relative Intensity')
plt.xlim([0,800])
plt.title('Polystyrene standard Size = 495 nm, CW = 200 µm')
plt.grid()
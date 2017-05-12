# -*- coding: utf-8 -*-
"""
Created on Sat Jan 28 18:22:34 2017

@author: Michael Orrill
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
rcParams['font.size'] = 14 
#%% read in data
co2 = np.loadtxt('CO2_data.txt',skiprows=1)
#%% plot data
plt.figure(figsize=(11,6))
plt.plot(co2[:,0]*1e-3,co2[:,1],color=(0,0.5,1),linewidth=3)
plt.gca().invert_xaxis()
plt.xlabel('Thousands of years before present (0 = 2001)')
plt.ylabel('Atmospheric $\mathregular{CO_2}$ levels (ppm)')
plt.xlim((np.max(co2[:,0]*1e-3),np.min(co2[:,0]*1e-3)-10))
plt.ylim((170,400))
plt.tick_params(axis='both',which='major',labelsize=14)
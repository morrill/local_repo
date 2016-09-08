# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#%%
import numpy as np
import matplotlib.pyplot as plt

#%% Reading in data from txt files
temp = np.loadtxt('water_ice_9_7_16.txt', delimiter='\t', skiprows=18, usecols=(0,2))
water_ice = temp*(1/60, 1) - (0, temp[0,1])
temp = np.loadtxt('water_no_ice_9_7_16.txt', delimiter='\t', skiprows=18, usecols=(0,3))
water_no_ice = temp*(1/60, 1) - (0, temp[0,1])

temp = np.loadtxt('EG_ice_9_7_16.txt', delimiter='\t', skiprows=18, usecols=(0,2))
EG_ice = temp*(1/60, 1) - (0, temp[0,1])
temp = np.loadtxt('EG_no_ice_9_7_16.txt', delimiter='\t', skiprows=20, usecols=(0,2))
EG_no_ice = temp*(1/60, 1) - (0, temp[0,1])

#%% Plotting data
fig = plt.figure(figsize=(20,15))
ax1 = fig.add_subplot(211)
ax1.plot(water_ice[:,0], water_ice[:,1], water_no_ice[:,0], water_no_ice[:,1], 'r-')
plt.xlabel('Time [min]')
plt.ylabel('Normalized Temabsoluter [$^{\circ}$C]')
plt.xlim([0,np.amax([np.amax(water_ice),np.amax(water_no_ice)])])
plt.legend(["ice bath", "no ice bath"])
plt.title('Water')

ax2 = fig.add_subplot(212)
ax2.plot(EG_ice[:,0], EG_ice[:,1], EG_no_ice[:,0], EG_no_ice[:,1], 'r-')
plt.xlabel('Time [min]')
plt.ylabel('Normalized Temabsoluter [$^{\circ}$C]')
plt.xlim([0,np.amax([np.amax(EG_ice),np.amax(EG_no_ice)])])
plt.legend(["ice bath", "no ice bath"])
plt.title('Ethylene Glycol')

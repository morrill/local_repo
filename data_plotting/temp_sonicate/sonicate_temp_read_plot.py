# -*- coding: utf-8 -*-
"""
This is for plotting temperature data from sonication treatments

This is a temporary script file.
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size':18})
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
fig = plt.figure(figsize=(20,22.5))
ax1 = fig.add_subplot(311)
ax1.plot(water_ice[:,0], water_ice[:,1], water_no_ice[:,0], water_no_ice[:,1], 'r-')
plt.xlabel('Time [min]')
plt.ylabel('Relative Temperature Change [$^{\circ}$C]')
plt.xlim([0,25])
plt.legend(["ice bath", "no ice bath"],loc=0)
plt.title('Water')

ax2 = fig.add_subplot(312)
ax2.plot(EG_ice[:,0], EG_ice[:,1], EG_no_ice[:,0], EG_no_ice[:,1], 'r-')
plt.xlabel('Time [min]')
plt.ylabel('Relative Temperature Change [$^{\circ}$C]')
plt.xlim([0,25])
plt.legend(["ice bath", "no ice bath"],loc=0)
plt.title('Ethylene Glycol')

ax2 = fig.add_subplot(313)
ax2.plot(Ag_EMD5730_ice[:,0], Ag_EMD5730_ice[:,1], 'b-', Ag_EMD5800_ice[:,0], Ag_EMD5800_ice[:,1], 'r-', Cu_ice[:,0], Cu_ice[:,1], 'g-')
plt.xlabel('Time [min]')
plt.ylabel('Relative Temperature Change [$^{\circ}$C]')
plt.xlim([0,25])
plt.legend(["Ag EMD5730", "Ag EMD5800", "Cu"],loc=0)
plt.title('Metallic inks in ice')

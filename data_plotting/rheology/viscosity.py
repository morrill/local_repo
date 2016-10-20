# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 19:17:58 2016

@author: Michael Orrill

This is for plotting results from rheometer viscosity measurements
"""
import numpy as np
import matplotlib.pyplot as plt
#%%
rt5 = np.loadtxt('RT5_standard.txt', skiprows=3, usecols=(0,1,2)) # RT5 Standard
eg11 = np.loadtxt('PureEG_1_1.txt', skiprows=3, usecols=(0,1,2)) # pure ethylene glycol 1.1
eg12 = np.loadtxt('PureEG_1_2.txt', skiprows=3, usecols=(0,1,2)) # pure ethylene glycol 1.2
eg21 = np.loadtxt('PureEG_2_1.txt', skiprows=3, usecols=(0,1,2)) # pure ethylene glycol 2.1
eg22 = np.loadtxt('PureEG_2_2.txt', skiprows=3, usecols=(0,1,2)) # pure ethylene glycol 2.2
C_son_11 = np.loadtxt('C_ink_sonicated1_1.txt', skiprows=3, usecols=(0,1,2)) # sonicated C ink 1.1
C_son_12 = np.loadtxt('C_ink_sonicated1_2.txt', skiprows=3, usecols=(0,1,2)) # sonicated C ink 1.2
C_son_21 = np.loadtxt('C_ink_sonicated2_1.txt', skiprows=3, usecols=(0,1,2)) # sonicated C ink 2.1
C_son_22 = np.loadtxt('C_ink_sonicated2_2.txt', skiprows=3, usecols=(0,1,2)) # sonicated C ink 2.2
C11 = np.loadtxt('C_ink1_1.txt', skiprows=3, usecols=(0,1,2)) # C ink 1.1
C12 = np.loadtxt('C_ink1_2.txt', skiprows=3, usecols=(0,1,2)) # C ink 1.2
C21 = np.loadtxt('C_ink2_1.txt', skiprows=3, usecols=(0,1,2)) # C ink 2.1
C22 = np.loadtxt('C_ink2_2.txt', skiprows=3, usecols=(0,1,2)) # C ink 2.2
#%% plot RT5 sample and calculate % error for average of semi-flat region 
rt5_ave = np.mean(rt5[np.min(np.where(rt5[:,1]>=2e3)):,2]*1e3)
rt5_data = np.array([[20,4.941],[23,4.682],[24,4.600],[25,4.520]]) # data from RT5 bottle
test_temp = 23.5 # testing temp
rt5_expect_visc = rt5_data[1,1] + (test_temp-rt5_data[1,0])*(rt5_data[2,1] - rt5_data[1,1])/(rt5_data[2,0] - rt5_data[1,0]) #interpolated expected viscosity
rt5_error = abs((rt5_ave - rt5_expect_visc)/rt5_expect_visc) # percent error

plt.figure(figsize=(10,5))
plt.plot(rt5[:,1],rt5[:,2]*1000, color=(0.8,0.5,0.1), label='Measured')
plt.axhline(rt5_expect_visc, label='Standard')
plt.xlabel('Shear rate [1/s]')
plt.ylabel('Viscosity [mPa.s]')
plt.legend()
plt.xlim([0,6000])
plt.ylim(0,30)
plt.title('RT5 standard, Error = %0.2f %%' % (rt5_error*100))
plt.grid()
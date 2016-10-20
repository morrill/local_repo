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
c_son_11 = np.loadtxt('C_ink_sonicated1_1.txt', skiprows=3, usecols=(0,1,2)) # sonicated C ink 1.1
c_son_12 = np.loadtxt('C_ink_sonicated1_2.txt', skiprows=3, usecols=(0,1,2)) # sonicated C ink 1.2
c_son_21 = np.loadtxt('C_ink_sonicated2_1.txt', skiprows=3, usecols=(0,1,2)) # sonicated C ink 2.1
c_son_22 = np.loadtxt('C_ink_sonicated2_2.txt', skiprows=3, usecols=(0,1,2)) # sonicated C ink 2.2
c11 = np.loadtxt('C_ink1_1.txt', skiprows=3, usecols=(0,1,2)) # C ink 1.1
c12 = np.loadtxt('C_ink1_2.txt', skiprows=3, usecols=(0,1,2)) # C ink 1.2
c21 = np.loadtxt('C_ink2_1.txt', skiprows=3, usecols=(0,1,2)) # C ink 2.1
c22 = np.loadtxt('C_ink2_2.txt', skiprows=3, usecols=(0,1,2)) # C ink 2.2
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
#%% plot EG sample and calculate % error for average of semi-flat region 
# get maximum error from nominal value
eg_expect_visc = 16.1 # data from several online sources at 25°C
eg_ave = np.zeros(4) # average valuse for each measurement
eg_ave[0] = np.mean(eg11[np.min(np.where(eg11[:,1]>=2e3)):,2]*1e3)
eg_ave[1] = np.mean(eg12[np.min(np.where(eg12[:,1]>=2e3)):,2]*1e3)
eg_ave[2] = np.mean(eg21[np.min(np.where(eg21[:,1]>=2e3)):,2]*1e3)
eg_ave[3] = np.mean(eg22[np.min(np.where(eg22[:,1]>=2e3)):,2]*1e3)

eg_error = np.zeros(4)
eg_error[0] = abs((eg_ave[0] - eg_expect_visc)/eg_expect_visc) # percent error from 1.1
eg_error[1] = abs((eg_ave[1] - eg_expect_visc)/eg_expect_visc) # percent error from 1.2
eg_error[2] = abs((eg_ave[2] - eg_expect_visc)/eg_expect_visc) # percent error from 2.1
eg_error[3] = abs((eg_ave[3] - eg_expect_visc)/eg_expect_visc) # percent error from 2.2
eg_max_error = eg_error.max()

plt.figure(figsize=(10,5))
plt.plot(eg11[:,1],eg11[:,2]*1000, color=(1,0,0), label='1.1') # pure ethylene glycol 1.1
plt.plot(eg12[:,1],eg12[:,2]*1000, color=(1,0.5,0), label='1.2') # pure ethylene glycol 1.2
plt.plot(eg21[:,1],eg21[:,2]*1000, color=(0,1,0), label='2.1') # pure ethylene glycol 2.1
plt.plot(eg22[:,1],eg22[:,2]*1000, color=(0,1,0.6), label='2.2') # pure ethylene glycol 2.2
plt.axhline(eg_expect_visc, label='Refernce at 25°C')
plt.xlabel('Shear rate [1/s]')
plt.ylabel('Viscosity [mPa.s]')
plt.legend()
plt.xlim([0,6000])
plt.ylim(10,25)
plt.title('Pure ethylene glycol (at 23.5 °C), Maximum Error = %0.2f %%' % (eg_max_error*100))
plt.grid()
#%% plot sonicated carbon ink sample 
plt.figure(figsize=(10,5))
plt.plot(c_son_11[:,1],c_son_11[:,2]*1000, color=(1,0,0), label='1.1') # sonicated C ink 1.1
plt.plot(c_son_12[:,1],c_son_12[:,2]*1000, color=(1,0.5,0), label='1.2') # sonicated C ink 1.2
plt.plot(c_son_21[:,1],c_son_21[:,2]*1000, color=(0,1,0), label='2.1') # sonicated C ink 2.1
plt.plot(c_son_22[:,1],c_son_22[:,2]*1000, color=(0,1,0.6), label='2.2') # sonicated C ink 2.2
plt.xlabel('Shear rate [1/s]')
plt.ylabel('Viscosity [mPa.s]')
plt.legend()
plt.xlim([0,6000])
plt.ylim(5,20)
plt.title('Sonicated Carbon ink (at 23.5 °C)')
plt.grid()
#%% plot un-sonicated carbon ink sample 
plt.figure(figsize=(10,5))
plt.plot(c11[:,1],c11[:,2]*1000, color=(1,0,0), label='1.1') # C ink 1.1
plt.plot(c12[:,1],c12[:,2]*1000, color=(1,0.5,0), label='1.2') # C ink 1.2
plt.plot(c21[:,1],c21[:,2]*1000, color=(0,1,0), label='2.1') # C ink 2.1
plt.plot(c22[:,1],c22[:,2]*1000, color=(0,1,0.6), label='2.2') # C ink 2.2
plt.xlabel('Shear rate [1/s]')
plt.ylabel('Viscosity [mPa.s]')
plt.legend()
plt.xlim([0,6000])
plt.ylim(10,25)
plt.title('Un-Sonicated Carbon ink (at 23.5 °C)')
plt.grid()
#%% plot un-sonicated vs sonicated carbon ink sample 
plt.figure(figsize=(10,5))
plt.plot(c11[:,1],c11[:,2]*1000, color=(1,0,0), label='Un-sonicated') # C ink 1.1
plt.plot(c12[:,1],c12[:,2]*1000, color=(1,0,0)) # C ink 1.2
plt.plot(c21[:,1],c21[:,2]*1000, color=(1,0,0)) # C ink 2.1
plt.plot(c22[:,1],c22[:,2]*1000, color=(1,0,0)) # C ink 2.2

plt.plot(c_son_11[:,1],c_son_11[:,2]*1000, color=(0,1,0), label='Sonicated') # sonicated C ink 1.1
plt.plot(c_son_12[:,1],c_son_12[:,2]*1000, color=(0,1,0)) # sonicated C ink 1.2
plt.plot(c_son_21[:,1],c_son_21[:,2]*1000, color=(0,1,0)) # sonicated C ink 2.1
plt.plot(c_son_22[:,1],c_son_22[:,2]*1000, color=(0,1,0)) # sonicated C ink 2.2

plt.xlabel('Shear rate [1/s]')
plt.ylabel('Viscosity [mPa.s]')
plt.legend()
plt.xlim([0,6000])
plt.ylim(5,25)
plt.title('Un-Sonicated vs sonicated Carbon ink (at 23.5 °C)')
plt.grid()

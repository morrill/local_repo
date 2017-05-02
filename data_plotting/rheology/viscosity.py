# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 19:17:58 2016

@author: Michael Orrill

This is for plotting results from rheometer viscosity measurements
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib import rcParams
rcParams['font.size'] = 14

#%% Read in data from 2-27-17
head_skip = np.arange(21,344,41)
for i, num in enumerate(head_skip):
    exec("cs%s = np.genfromtxt('c4.1_sonicated.txt', skip_header=%d, max_rows=36, usecols=(0,1,2))" % (i,num))
    print(i)
#%% Plot data from 2-27-17
yellow=np.linspace(0,0.9,8)
fig = plt.figure(figsize=(11,6))
ax = fig.add_subplot(111)
for i, num in enumerate(range(0,8)):
    exec("ax.semilogx(cs%s[:,1],cs%s[:,2]*1000, color=(1,yellow[%d],0))" % (num,num,i))
plt.xlabel('Shear rate [$s^{-1}$]')
plt.ylabel('Viscosity [$mPa*s$]')
plt.xlim([10,2500])
plt.ylim(6,18)
ax.xaxis.set_major_formatter(ScalarFormatter())
plt.title('c4.1, 8 measurements, solvent trap, sonication')
plt.tick_params(labelright=True)
#%% Plot averages from 2-27-17
ave = np.zeros(8)
for i in range(0,8):
    exec("ave[%d] = np.mean(cs%s[:,2]*1000)" % (i,i))
plt.figure(figsize=(11,6))
plt.plot(np.arange(1,9),ave,color=(1,0.2,0),lw=3)
plt.xlabel('Measurement number')
plt.ylabel('Viscosity [$mPa*s$]')
plt.ylim(6,18)
plt.xlim([1,10])
plt.title('Average viscosity of each measurement, c4.1 w/ solvent trap, sonication')
#%% Statistics from 2-27-17
sd = np.zeros(8)
for i in range(0,8):
    exec("sd[%d] = np.std(cs%s[:,2]*1000)" % (i,i))

#%% Read in data from 2-24-17
head_skip = np.arange(21,385,41)
for i, num in enumerate(head_skip):
    exec("c%s = np.genfromtxt('C4.1.txt', skip_header=%d, max_rows=36, usecols=(0,1,2))" % (i,num))
#%% Plot data from 2-24-17
yellow=np.linspace(0,0.9,9)
fig = plt.figure(figsize=(11,6))
ax = fig.add_subplot(111)
for i, num in enumerate(range(0,9)):
    exec("ax.semilogx(c%s[:,1],c%s[:,2]*1000, color=(1,yellow[%d],0))" % (num,num,i))
plt.xlabel('Shear rate [$s^{-1}$]')
plt.ylabel('Viscosity [$mPa*s$]')
plt.xlim([10,2500])
plt.ylim(6,18)
ax.xaxis.set_major_formatter(ScalarFormatter())
plt.title('C4.1, 9 measurements, solvent trap, no sonication')
plt.tick_params(labelright=True)
#%% Plot averages from 2-24-17
ave = np.zeros(9)
for i in range(0,9):
    exec("ave[%d] = np.mean(c%s[:,2]*1000)" % (i,i))
plt.figure(figsize=(11,6))
plt.plot(np.arange(1,10),ave,color=(1,0.2,0),lw=3)
plt.xlabel('Measurement number')
plt.ylabel('Viscosity [$mPa*s$]')
plt.ylim(6,18)
plt.xlim([1,10])
plt.title('Average viscosity of each measurement, C4.1 w/ solvent trap, no sonication')



#%% plot all from 2-24-17 and 2-27-17 and 2-17-17
a1 = np.linspace(0,0.9,9)[::-1]
fig = plt.figure(figsize=(11,6))
ax = fig.add_subplot(111)
for i, num in enumerate(range(0,1)):
    exec("ax.semilogx(c%s[:,1],c%s[:,2]*1000, color='k', alpha=a1[i], label='C4.1 Not sonicated' if i == 0 else '')" % (num,num))
for i, num in enumerate(range(0,1)):  
    exec("ax.semilogx(cs%s[:,1],cs%s[:,2]*1000, color='b', alpha=a1[i], label='C4.1 Sonicated'if i == 0 else '')" % (num,num))
for i, num in enumerate(range(0,1)):  
    exec("ax.semilogx(eg%s[:,1],eg%s[:,2]*1000, color='g', alpha=a1[i], label='Neat ethylene glycol'if i == 0 else '')" % (num,num))
plt.xlabel(r'Shear rate ($\mathrm{s}^{-1}$)')
plt.ylabel(r'Viscosity ($\mathrm{mPa*s}$)')
plt.xlim([10,2500])
plt.ylim(6,18)
ax.xaxis.set_major_formatter(ScalarFormatter())
#plt.title('C4.1, 9 measurements, solvent trap, no sonication')
plt.tick_params(labelright=True)
plt.legend(loc=4)

#%% Plot averages from 2-24-17 and 2-27-17 and 2-17-17
avec = np.zeros(9)
avecs = np.zeros(8)
for i in range(0,9):
    exec("avec[%d] = np.mean(c%s[:,2]*1000)" % (i,i))
for i in range(0,8):    
    exec("avecs[%d] = np.mean(cs%s[:,2]*1000)" % (i,i))
plt.figure(figsize=(11,6))
plt.plot(np.arange(1,9),avec[:8],color='k',lw=3,label='C4.1 Not sonicated')
plt.plot(np.arange(1,9),avecs,color='b',lw=3,label='C4.1 Sonicated')
plt.plot(np.arange(1,9),ave[:8],color='g',lw=3,label='Neat ethylene glycol')
plt.xlabel('Measurement number')
plt.ylabel('Viscosity [$mPa*s$]')
plt.ylim(6,18)
plt.xlim([1,9])
plt.title('Average viscosity of each measurement, C4.1 w/ solvent trap')
plt.legend(loc=4)

#%% Read in data from 2-17-17
head_skip = np.arange(27,472,30)
for i, num in enumerate(head_skip):
    exec("eg%s = np.genfromtxt('EG.txt', skip_header=%d, max_rows=25, usecols=(0,1,2))" % (i,num))
#%% Plot data from 2-17-17
yellow=np.linspace(0,0.9,14)
fig = plt.figure(figsize=(11,6))
ax = fig.add_subplot(111)
for i, num in enumerate(range(19,35)):
    exec("ax.semilogx(eg%s[:,1],eg%s[:,2]*1000, color=(1,yellow[%d],0))" % (num,num,i))
plt.xlabel('Shear rate [$s^{-1}$]')
plt.ylabel('Viscosity [$mPa*s$]')
plt.xlim([10,12000])
plt.ylim(6,18)
ax.xaxis.set_major_formatter(ScalarFormatter())
plt.title('Ethylene glycol, last 16 measurements')
plt.tick_params(labelright=True)
#%% Plot averages from 2-17-17
ave = np.zeros(14)
for i in range(0,14):
    exec("ave[%d] = np.mean(eg%s[:,2]*1000)" % (i,i))
plt.figure(figsize=(11,6))
plt.plot(np.arange(1,15),ave,color=(1,0.2,0),lw=3)
plt.xlabel('Measurement number')
plt.ylabel('Viscosity [$mPa*s$]')
plt.ylim(6,18)
plt.xlim([1,35])
plt.title('Average viscosity of each measurement, EG w/ solvent trap')





#%% Read in data from 11-30-16
head_skip = np.arange(46,1305,37)
for i, num in enumerate(head_skip):
    exec("eg%s = np.genfromtxt('EG.txt', skip_header=%d, max_rows=31, usecols=(0,1,2))" % (i,num))
#%% Plot data from 11-30-16
yellow=np.linspace(0,0.9,16)
fig = plt.figure(figsize=(11,6))
ax = fig.add_subplot(111)
for i, num in enumerate(range(0,14)):
    exec("ax.semilogx(eg%s[:,1],eg%s[:,2]*1000, color=(1,yellow[%d],0))" % (num,num,i))
plt.xlabel('Shear rate [$s^{-1}$]')
plt.ylabel('Viscosity [$mPa*s$]')
plt.xlim([10,2500])
plt.ylim(13.5,14.5)
ax.xaxis.set_major_formatter(ScalarFormatter())
plt.title('Ethylene glycol, 14 measurements w/ solvent trap')
plt.tick_params(labelright=True)
#%% Plot averages from 11-30-16
ave = np.zeros(35)
for i in range(0,35):
    exec("ave[%d] = np.mean(eg%s[5:,2]*1000)" % (i,i))
plt.figure(figsize=(11,6))
plt.plot(np.arange(1,36),ave,color=(1,0.2,0),lw=3)
plt.xlabel('Measurement number')
plt.ylabel('Viscosity [$mPa*s$]')
plt.ylim(6,18)
plt.xlim([1,35])
plt.title('Average viscosity of each measurement')
#%% Read Data from 11-15-16
#rt5 = np.loadtxt('RT5_standard.txt', skiprows=3, usecols=(0,1,2)) # RT5 Standard
eg11 = np.loadtxt('eg11.txt', skiprows=2, usecols=(0,1,2)) # pure ethylene glycol 1.1
eg12 = np.loadtxt('eg12.txt', skiprows=2, usecols=(0,1,2)) # pure ethylene glycol 1.2
eg13 = np.loadtxt('eg13.txt', skiprows=2, usecols=(0,1,2)) # pure ethylene glycol 1.3
eg14 = np.loadtxt('eg14.txt', skiprows=2, usecols=(0,1,2)) # pure ethylene glycol 1.4
eg15 = np.loadtxt('eg15.txt', skiprows=2, usecols=(0,1,2)) # pure ethylene glycol 1.5
eg21 = np.loadtxt('eg21.txt', skiprows=2, usecols=(0,1,2)) # pure ethylene glycol 2.1
eg22 = np.loadtxt('eg22.txt', skiprows=2, usecols=(0,1,2)) # pure ethylene glycol 2.2
eg23 = np.loadtxt('eg23.txt', skiprows=2, usecols=(0,1,2)) # pure ethylene glycol 2.3
eg24 = np.loadtxt('eg24.txt', skiprows=2, usecols=(0,1,2)) # pure ethylene glycol 2.4
eg25 = np.loadtxt('eg25.txt', skiprows=2, usecols=(0,1,2)) # pure ethylene glycol 2.5
eg26 = np.loadtxt('eg26.txt', skiprows=2, usecols=(0,1,2)) # pure ethylene glycol 2.6
eg27 = np.loadtxt('eg27.txt', skiprows=2, usecols=(0,1,2)) # pure ethylene glycol 2.7
eg28 = np.loadtxt('eg28.txt', skiprows=2, usecols=(0,1,2)) # pure ethylene glycol 2.8
eg29 = np.loadtxt('eg29.txt', skiprows=2, usecols=(0,1,2)) # pure ethylene glycol 2.9
#%% plot Data from 11-15-16
yellow=np.linspace(0,0.9,9)
fig = plt.figure(figsize=(11,6))
ax = fig.add_subplot(111)
#ax.semilogx(eg11[:,1],eg11[:,2]*1000, color=(1,yellow[0],0), label='2.1')
#ax.semilogx(eg12[:,1],eg12[:,2]*1000, color=(1,yellow[1],0), label='2.2')
#ax.semilogx(eg13[:,1],eg13[:,2]*1000, color=(1,yellow[2],0), label='2.3')
#ax.semilogx(eg14[:,1],eg14[:,2]*1000, color=(1,yellow[3],0), label='2.4')
#ax.semilogx(eg15[:,1],eg15[:,2]*1000, color=(1,yellow[4],0), label='2.5')

ax.semilogx(eg21[:,1],eg21[:,2]*1000, color=(1,yellow[0],0), label='2.1')
ax.semilogx(eg22[:,1],eg22[:,2]*1000, color=(1,yellow[1],0), label='2.2')
ax.semilogx(eg23[:,1],eg23[:,2]*1000, color=(1,yellow[2],0), label='2.3')
ax.semilogx(eg24[:,1],eg24[:,2]*1000, color=(1,yellow[3],0), label='2.4')
ax.semilogx(eg25[:,1],eg25[:,2]*1000, color=(1,yellow[4],0), label='2.5')
ax.semilogx(eg26[:,1],eg26[:,2]*1000, color=(1,yellow[5],0), label='2.6')
ax.semilogx(eg27[:,1],eg27[:,2]*1000, color=(1,yellow[6],0), label='2.7')
ax.semilogx(eg28[:,1],eg28[:,2]*1000, color=(1,yellow[7],0), label='2.8')
ax.semilogx(eg29[:,1],eg29[:,2]*1000, color=(1,yellow[8],0), label='2.9')

plt.axhline(16.1, label='Spec') # horizontal line for the expected measurement
plt.xlabel('Shear rate [$s^{-1}$]')
plt.ylabel('Viscosity [$mPa*s$]')
plt.legend(loc='upper right',ncol=2, fontsize=14)
plt.xlim([1,12000])
ax.xaxis.set_major_formatter(ScalarFormatter())
plt.ylim(10,20)
plt.title('Ethylene Glycol at 23.0 °C (boiled w/ water near system)')
plt.tick_params(labelright=True)
#%% Read Data from 11-14-16
rt11 = np.loadtxt('RT511.txt', skiprows=3, usecols=(0,1,2)) # RT5 Standard
rt12 = np.loadtxt('RT512.txt', skiprows=3, usecols=(0,1,2)) # RT5 Standard
rt21 = np.loadtxt('RT521.txt', skiprows=3, usecols=(0,1,2)) # RT5 Standard
rt22 = np.loadtxt('RT522.txt', skiprows=3, usecols=(0,1,2)) # RT5 Standard
rt23 = np.loadtxt('RT523.txt', skiprows=3, usecols=(0,1,2)) # RT5 Standard
rt24 = np.loadtxt('RT524.txt', skiprows=3, usecols=(0,1,2)) # RT5 Standard
rt25 = np.loadtxt('RT525.txt', skiprows=3, usecols=(0,1,2)) # RT5 Standard
#%% Plot data from 11-14-16
yellow=np.linspace(0,0.9,5)
fig = plt.figure(figsize=(11,6))
ax = fig.add_subplot(111)

#ax.semilogx(rt11[:,1],rt11[:,2]*1000, color=(1,yellow[0],0), label='1.1')
#ax.semilogx(rt12[:,1],rt12[:,2]*1000, color=(1,yellow[3],0), label='1.2')

ax.semilogx(rt21[:,1],rt21[:,2]*1000, color=(1,yellow[0],0), label='2.1')
ax.semilogx(rt22[:,1],rt22[:,2]*1000, color=(1,yellow[1],0), label='2.2')
ax.semilogx(rt23[:,1],rt23[:,2]*1000, color=(1,yellow[2],0), label='2.3')
ax.semilogx(rt24[:,1],rt24[:,2]*1000, color=(1,yellow[3],0), label='2.2')
ax.semilogx(rt25[:,1],rt25[:,2]*1000, color=(1,yellow[4],0), label='2.3')

plt.axhline(4.682, label='Spec') # horizontal line for the expected measurement
plt.xlabel('Shear rate [$s^{-1}$]')
plt.ylabel('Viscosity [$mPa*s$]')
plt.legend(loc='upper right',ncol=1, fontsize=14)
plt.xlim([1,7000])
ax.xaxis.set_major_formatter(ScalarFormatter())
plt.ylim(1,10)
plt.title('RT5 Standard at 23.0 °C')
plt.tick_params(labelright=True)
#%% plot sonicated carbon ink sample 
yellow=np.linspace(0,0.9,9)

fig = plt.figure(figsize=(11,6))
ax = fig.add_subplot(111)
#ax.semilogx(eg11[:,1],eg11[:,2]*1000, color=(1,yellow[0],0), label='2.1')
#ax.semilogx(eg12[:,1],eg12[:,2]*1000, color=(1,yellow[1],0), label='2.2')
#ax.semilogx(eg13[:,1],eg13[:,2]*1000, color=(1,yellow[2],0), label='2.3')
#ax.semilogx(eg14[:,1],eg14[:,2]*1000, color=(1,yellow[3],0), label='2.4')
#ax.semilogx(eg15[:,1],eg15[:,2]*1000, color=(1,yellow[4],0), label='2.5')




ax.semilogx(eg21[:,1],eg21[:,2]*1000, color=(1,yellow[0],0), label='2.1')
ax.semilogx(eg22[:,1],eg22[:,2]*1000, color=(1,yellow[1],0), label='2.2')
ax.semilogx(eg23[:,1],eg23[:,2]*1000, color=(1,yellow[2],0), label='2.3')
ax.semilogx(eg24[:,1],eg24[:,2]*1000, color=(1,yellow[3],0), label='2.4')
ax.semilogx(eg25[:,1],eg25[:,2]*1000, color=(1,yellow[4],0), label='2.5')
ax.semilogx(eg26[:,1],eg26[:,2]*1000, color=(1,yellow[5],0), label='2.6')
ax.semilogx(eg27[:,1],eg27[:,2]*1000, color=(1,yellow[6],0), label='2.7')
ax.semilogx(eg28[:,1],eg28[:,2]*1000, color=(1,yellow[7],0), label='2.8')
ax.semilogx(eg29[:,1],eg29[:,2]*1000, color=(1,yellow[8],0), label='2.9')



plt.axhline(16.1, label='Spec')

#plt.plot(st11[:,1],st11[:,2]*1000, color=(1,0,0), label='1.1') 
#plt.plot(st12[:,1],st12[:,2]*1000, color=(1,0.2,0), label='1.2') 
#plt.plot(eg13[:,1],eg13[:,2]*1000, color=(1,0.4,0), label='1.3')
#plt.plot(eg14[:,1],eg14[:,2]*1000, color=(1,0.6,0), label='1.4')
#plt.plot(eg15[:,1],eg15[:,2]*1000, color=(1,0.8,0), label='1.5')
plt.xlabel('Shear rate [$s^{-1}$]')
plt.ylabel('Viscosity [$mPa*s$]')
plt.legend(loc='upper right',ncol=2, fontsize=14)
plt.xlim([1,12000])
ax.xaxis.set_major_formatter(ScalarFormatter())
plt.ylim(10,20)
plt.title('Ethylene Glycol at 23.0 °C (boiled w/ water near system)')
plt.tick_params(labelright=True)
#plt.xscale('log')

#plt.figure(figsize=(10,5))
#plt.plot(eg21[:,1],eg21[:,2]*1000, color=(1,0,0), label='2.1') 
#plt.plot(eg22[:,1],eg22[:,2]*1000, color=(1,0.2,0), label='2.2') 
#plt.plot(eg23[:,1],eg23[:,2]*1000, color=(1,0.4,0), label='2.3')
#plt.plot(eg24[:,1],eg24[:,2]*1000, color=(1,0.6,0), label='2.4')
#plt.plot(eg25[:,1],eg25[:,2]*1000, color=(1,0.8,0), label='2.5')
##plt.plot(eg31[:,1],eg31[:,2]*1000, color=(0,0,1), label='~2.6')
#plt.xlabel('Shear rate [1/s]')
#plt.ylabel('Viscosity [mPa.s]')
#plt.legend()
#plt.xlim([1,10000])
#plt.ylim(7,20)
#plt.title('Pure ethylene glycol (at 23.0 °C)')
##plt.xscale('log')

#plt.figure(figsize=(10,5))
#plt.plot(eg11[:,1],eg11[:,2]*1000, color=(1,0,0), label='1.1') 
#plt.plot(eg12[:,1],eg12[:,2]*1000, color=(1,0,0.8), label='1.2') 
#plt.plot(eg13[:,1],eg13[:,2]*1000, color=(0,1,0), label='1.3')
#plt.plot(eg14[:,1],eg14[:,2]*1000, color=(0,0,1), label='1.4')
#plt.plot(eg15[:,1],eg15[:,2]*1000, color=(0.5,0.5,0.5), label='1.5')
#plt.plot(eg21[:,1],eg21[:,2]*1000, color=(1,0,0), label='2.1') 
#plt.plot(eg22[:,1],eg22[:,2]*1000, color=(1,0,0.8), label='2.2') 
#plt.plot(eg23[:,1],eg23[:,2]*1000, color=(0,1,0), label='2.3')
#plt.plot(eg24[:,1],eg24[:,2]*1000, color=(0,0,1), label='2.4')
#plt.plot(eg25[:,1],eg25[:,2]*1000, color=(0.5,0.5,0.5), label='2.5')
#plt.xlabel('Shear rate [1/s]')
#plt.ylabel('Viscosity [mPa.s]')
#plt.legend(loc='upper right', ncol=2, fontsize=14)
#plt.xlim([1,10000])
#plt.ylim(9,17)
#plt.title('Pure ethylene glycol (at 23.0 °C)')
##plt.xscale('log')

#plt.figure(figsize=(10,5))
#plt.plot(eg15[:,1],eg15[:,2]*1000, color=(1,0,0), label='1.5')
#plt.plot(eg25[:,1],eg25[:,2]*1000, color=(0,1,0), label='2.5')
#plt.plot(eg31[:,1],eg31[:,2]*1000, color=(0,0,1), label='2.5')
#plt.xlabel('Shear rate [1/s]')
#plt.ylabel('Viscosity [mPa.s]')
#plt.legend()
#plt.xlim([1,10000])
#plt.ylim(7,13)
#plt.title('Pure ethylene glycol (at 23.0 °C)')
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
#%% read in data from 10-19-16
for samp in range(1,3):
    for meas in range(1,3):
        exec("c%s%s = np.loadtxt('c%s%s.txt',skiprows=3, usecols=(0,1,2))" % (samp,meas,samp,meas))
        exec("cs%s%s = np.loadtxt('cs%s%s.txt',skiprows=3, usecols=(0,1,2))" % (samp,meas,samp,meas))
#%% plot un-sonicated vs sonicated carbon ink sample 
plt.figure(figsize=(20,10))
plt.plot(c11[:,1],c11[:,2]*1000, color='b', linewidth=3, label='Un-sonicated') # C ink 1.1
plt.plot(c12[:,1],c12[:,2]*1000, color='b', linewidth=3) # C ink 1.2
plt.plot(c21[:,1],c21[:,2]*1000, color='b', linewidth=3) # C ink 2.1
plt.plot(c22[:,1],c22[:,2]*1000, color='b', linewidth=3) # C ink 2.2

plt.plot(cs11[:,1],cs11[:,2]*1000, color='k', linewidth=3, label='Sonicated') # sonicated C ink 1.1
plt.plot(cs12[:,1],cs12[:,2]*1000, color='k', linewidth=3) # sonicated C ink 1.2
plt.plot(cs21[:,1],cs21[:,2]*1000, color='k', linewidth=3) # sonicated C ink 2.1
plt.plot(cs22[:,1],cs22[:,2]*1000, color='k', linewidth=3) # sonicated C ink 2.2

plt.xlabel('Shear rate $\mathregular{(s^{-1})}$')
plt.ylabel('Viscosity (mPa s)')
plt.legend(fontsize=30)
plt.xlim([0,6000])
plt.ylim(6,25)
#plt.title('Un-Sonicated vs sonicated Carbon ink (at 23.5 °C)')
plt.grid()

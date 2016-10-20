# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 15:21:27 2016

@author: Michael Orrill
"""
#%%
import numpy as np
import matplotlib.pyplot as plt

<<<<<<< HEAD
#%% read in polystyrene standards
st_auto_off_200 = np.loadtxt('standard1_auto_channel_width_off.asc', delimiter="\t", skiprows=89)
st_auto_on1_200 = np.loadtxt('standard1_auto_channel_width_on.asc', delimiter="\t", skiprows=89)
st_auto_on2_200 = np.loadtxt('standard1_auto_channel_width_on2.asc', delimiter="\t", skiprows=89)
st1_CW_50um = np.loadtxt('standard1_CW_50µm.asc', delimiter="\t", skiprows=89)
st1_CW_70um = np.loadtxt('standard1_CW_70µm.asc', delimiter="\t", skiprows=89)
st1_CW_300um = np.loadtxt('standard1_CW_300µm.asc', delimiter="\t", skiprows=89) #*** channel width saved as 200 µm in file...
st1_CW_400um = np.loadtxt('standard1_CW_400µm.asc', delimiter="\t", skiprows=89)

#%% read in carbon ink from 9/1/16
c41 = np.loadtxt('C4.1_30ul_9_1_16.asc', delimiter="\t", skiprows=89)
c42 = np.loadtxt('C4.1_30ul_9_1_16.asc', delimiter="\t", skiprows=89)
c43 = np.loadtxt('C4.1_30ul_9_1_16.asc', delimiter="\t", skiprows=89)
=======
#%% read in all polystyrene standard measuerments
st0 = np.loadtxt('standard0.asc', delimiter="\t", skiprows=89)
st1 = np.loadtxt('standard1.asc', delimiter="\t", skiprows=89)
st2 = np.loadtxt('standard2.asc', delimiter="\t", skiprows=89)
st3 = np.loadtxt('standard3_9_6_16.asc', delimiter="\t", skiprows=89)
st4 = np.loadtxt('standard4_9_6_16.asc', delimiter="\t", skiprows=89)
st5 = np.loadtxt('standard5_9_6_16_new_water.asc', delimiter="\t", skiprows=89)
st6 = np.loadtxt('standard6_9_6_16_new_water.asc', delimiter="\t", skiprows=89)

#%% read in carbon ink from 9/1/16
c41 = np.loadtxt('C4.1_30ul_9_1_16.asc', delimiter="\t", skiprows=89)
c42 = np.loadtxt('C4.2_30ul_9_1_16.asc', delimiter="\t", skiprows=89)
c43 = np.loadtxt('C4.3_30ul_9_1_16.asc', delimiter="\t", skiprows=89)

#%% read in carbon ink from 9/6/16
c44 = np.loadtxt('C4.1_sonicated_9_6_16.asc', delimiter="\t", skiprows=89)

>>>>>>> 4fe93bbe4c4687abc7e0c741806b454fb2222e80
#%% Plot standards with a bar plot
plt.figure(figsize=(15,5))
plt.bar(st0[:,0],st0[:,1], width=15, color=(0.6,0.4,0.7,0.2))
plt.bar(st1[:,0],st1[:,1], width=15, color=(0,0.1,0.3,0.2))
plt.bar(st2[:,0],st2[:,1], width=15, color=(0,0.2,0.5,0.2))
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
plt.figure(figsize=(15,5))
<<<<<<< HEAD
plt.bar(st_auto_on2_200[:,0],st_auto_on2_200[:,1], width=15, color=(0.8,0.5,1,0.2))
plt.xlabel('Particle Size [nm]')
plt.ylabel('Relative Intensity')
plt.xlim([0,800])
plt.title('Polystyrene standard Size = 495 nm, CW = 200 µm')
=======
plt.subplot
plt.bar(c41[:,0],c41[:,1], width=5, color=(0,0.6,0.2,0.2))
plt.bar(c42[:,0],c42[:,1], width=5, color=(0,1,1,0.2))
plt.bar(c43[:,0],c43[:,1], width=5, color=(0,0.1,0.5,0.2))
# sonicated carbon sample from 9/6/16
plt.bar(c44[:,0],c44[:,1], width=5, color=(1,0.6,0.2,0.5))

plt.xlabel('Particle Size [nm]')
plt.ylabel('Relative Intensity')
plt.xlim([0,250])
plt.title('Carbon sample 4, 9/1/16 with sonicated sample on 9/6/16')
>>>>>>> 4fe93bbe4c4687abc7e0c741806b454fb2222e80
plt.grid()

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
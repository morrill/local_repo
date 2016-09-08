# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 11:54:17 2016

@author: Michael Orrill

TGA plotting
"""
#%%
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#%% Read in data and plot

d1 = np.loadtxt("Ag1_8_31_16_%vtemp.txt", delimiter="\t", skiprows=3)
d2 = np.loadtxt("Ag2_8_31_16_%vtemp.txt", delimiter="\t", skiprows=3)
d3 = np.loadtxt("Ag1_8_31_16_%vtime.txt", delimiter="\t", skiprows=3)
d4 = np.loadtxt("Ag2_8_31_16_%vtime.txt", delimiter="\t", skiprows=3)
  
fig =plt.figure(figsize=(16,10))
ax1 = fig.add_subplot(2,2,1)
ax2 = fig.add_subplot(2,2,2)
ax3 = fig.add_subplot(2,2,3)
ax4 = fig.add_subplot(2,2,4)
    
ax1.plot(d1[:,0],d1[:,1], sns.xkcd_rgb["strawberry"])
ax1.set_xlabel('Temperature [$^\circ$C]')
ax1.set_ylabel('Mass percentage [%]')
ax1.set_title('EMD 5730, Ag ink 1: solvent based')

ax2.plot(d2[:,0],d2[:,1], sns.xkcd_rgb["strawberry"])
ax2.set_xlabel('Temperature [$^\circ$C]')
ax2.set_ylabel('Mass percentage [%]')
ax2.set_title('EMD 5800, Ag ink 2: oil based')

ax3.plot(d3[:,0],d3[:,1], sns.xkcd_rgb["strawberry"])
ax3.set_xlabel('Time [min]')
ax3.set_ylabel('Mass percentage [%]')
ax3.set_title('EMD 5730, Ag ink 1: solvent based')

ax4.plot(d4[:,0],d4[:,1], sns.xkcd_rgb["strawberry"])
ax4.set_xlabel('Time [min]')
ax4.set_ylabel('Mass percentage [%]')
ax4.set_title('EMD 5800, Ag ink 2: oil based')

plt.tight_layout()

    
    
    

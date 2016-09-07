# -*- coding: utf-8 -*-
"""
Created on Thu May 05 11:52:57 2016

@author: Michael Orrill
The arrays needed to plot Max power density and efficiency vs rteo are saved in the plot data folder
"""
#%% Libraries
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import math as m
#%% Material properties:
# temperature
Th = 372.0        # hot stream tempertaure [K]
Tc = 294.0         # cold stream temperature [K]
# geometry
rwi = 3e-3        # inner pipe radius [m]
rwo = 4e-3         # outer pipe
tc = 100e-6         #contact width
ts = 100e-6          # spacer thickness for both radial and axial spacers [m]
tw = 10e-3       # pipe thickness [m]
rtei = rwo + ts #inner TE radius [m]
# material properties
kw = 385.0                  # pipe thermal conductivity [W/m.K]
ks = 0.1                # spacer/insulator thermal conductivity [W/m.K]
Rc = 1e-9                    # contact resistance [Ohm.m^2]
# thermal resistances for first two nodes:
Rw = m.log(rwo/rwi)/(2.0*m.pi*tw*kw)   # inner pipe thermal resistance
Rs = m.log(rtei/rwo)/(2.0*m.pi*tw*ks)
Rwpluss = Rw + Rs
#%% nanoplatlet material properties
kn = 0.51                 # n-type material thermal conductivity [W/m.K]
kp = 0.46                 # p-type material thermal conductivity [W/m.K]
Sn = -193e-6                 # n-type Seebeck coefficient [V/K] 
Sp = 225e-6                  # p-type Seebeck coefficient [V/K] 
sigma_n = 39000            # n-type electrical conductivity [S/m]
sigma_p = 35000            # p-type electrical conductivity [S/m]
#%% nanomaterial material properties
kn = 1.15                 # n-type material thermal conductivity [W/m.K]
kp = 1.05                 # p-type material thermal conductivity [W/m.K]
Sn = -200e-6                 # n-type Seebeck coefficient [V/K] 
Sp = 195e-6                  # p-type Seebeck coefficient [V/K] 
sigma_n = 87500            # n-type electrical conductivity [S/m]
sigma_p = 110000            # p-type electrical conductivity [S/m]
#%% bulk material properties
kn = 1.3                 # n-type material thermal conductivity [W/m.K]
kp = 1.35                 # p-type material thermal conductivity [W/m.K]
Sn = -150e-6                 # n-type Seebeck coefficient [V/K] 
Sp = 225e-6                  # p-type Seebeck coefficient [V/K] 
sigma_n = 86957            # n-type electrical conductivity [S/m]
sigma_p = 90000            # p-type electrical conductivity [S/m]
#%% plots for sliced data to give more accurate curves of of power density vs. leg length by finding
#mm values and tn values that maximizes power density for each value of rteo

# need to run this and raname pmax_rteo to start with "material name"_pmax_rteo etc.
mm_maxp_rteo = np.zeros(len(rteo))
tn_maxp_rteo = np.zeros(len(rteo))
rteo_maxp = np.zeros(len(rteo))
pmax_rteo_tn = np.zeros((len(rteo), len(tn)))
pmax_rteo = np.zeros(len(rteo))
tn_maxp = np.zeros(len(rteo))
hc_maxp = np.zeros(len(rteo))

for i in range(len(rteo)):
    [mm_maxp_rteo[i],tn_maxp_rteo[i]] = np.unravel_index(p[:,i,:].argmax(), p[:,i,:].shape) # gets rteo and tn for max p in each mm slice
    rteo_maxp[i] = p[:,i,:].max() # gets max p value from each mm value

for i in range(len(rteo)):
    pmax_rteo_tn[i, :] = p[int(round(mm_maxp_rteo[i])),i,:]
    pmax_rteo[i] = p[int(round(mm_maxp_rteo[i])),i, int(round(tn_maxp_rteo[i]))]
    tn_maxp[i] = tn[int(tn_maxp_rteo[i])]
    hc_maxp[i] = outconv[int(round(mm_maxp_rteo[i])),i, int(round(tn_maxp_rteo[i]))]
    
maxp2 = rteo_maxp.max()        # overall max power
maxprteo2 = rteo[int(rteo_maxp.argmax())]
maxpmm2 = mm[int(mm_maxp_rteo[rteo_maxp.argmax()])]
maxptn2 = tn[int(tn_maxp_rteo[rteo_maxp.argmax()])]
#%% Slices throuogh eta arrays to get values at maximum efficiency

# change the name of bulk_eta to whatever 
mm_maxeta_rteo = np.zeros(len(rteo))
tn_maxeta_rteo = np.zeros(len(rteo))
rteo_maxeta = np.zeros(len(rteo))
maxeta_rteo_tn = np.zeros((len(rteo), len(tn)))
maxeta_rteo = np.zeros(len(rteo))

for i in range(len(rteo)):
    [mm_maxeta_rteo[i],tn_maxeta_rteo[i]] = np.unravel_index(bulk_eta[:,i,:].argmax(), bulk_eta[:,i,:].shape) # gets rteo and tn for maxeta in each mm slice
    rteo_maxeta[i] = bulk_eta[:,i,:].max() # gets max p value from each mm value

for i in range(len(rteo)):
    maxeta_rteo_tn[i, :] = bulk_eta[int(round(mm_maxeta_rteo[i])),i,:]
    maxeta_rteo[i] = bulk_eta[int(round(mm_maxeta_rteo[i])),i, int(round(tn_maxeta_rteo[i]))]

maxeta2 = rteo_maxeta.max()        # overall max power
maxetarteo2 = rteo[int(rteo_maxeta.argmax())]
maxetamm2 = mm[int(mm_maxeta_rteo[rteo_maxeta.argmax()])]
maxetatn2 = tn[int(tn_maxeta_rteo[rteo_maxeta.argmax()])]
#%% get efficiency
nano_Qh = (nano_Tw-nano_T1)/Rwpluss
bulk_Qh = (bulk_Tw-bulk_T1)/Rwpluss
nanoplt_Qh = (nanoplt_Tw-nanoplt_T1)/Rwpluss

nano_eta = (nano_p*2.0*m.pi*rtei*tw)/nano_Qh
bulk_eta = (bulk_p*2.0*m.pi*rtei*tw)/bulk_Qh
nanoplt_eta = (nanoplt_p*2.0*m.pi*rtei*tw)/nanoplt_Qh
#%% Get thermal and electrical resistances

Rte = (((2.0*m.pi)/(np.log(rteo/rtei)))*(kn*tn_maxp+ks*ts+kp*tp_maxp))**(-1.0) # TE/spacer total thermal resistance
Relec = ((np.log(rteo/rtei))/(2.0*m.pi*sigma_n*tn_maxp)) + ((np.log(rteo/rtei))/(2.0*m.pi*sigma_p*tp_maxp)) + (Rc/((2*m.pi*((rtei+tc)**2.0 - (rtei)**2.0)) + (2*m.pi*(((rteo)**2.0) - (rteo-tc)**2.0)))) # n- and p-type total electrical resistance plus contact resitance




#%% plot things
plt.figure(1)
#mat_name = ['Nanomaterial','Bulk','Nanoplatlets'] 
mat_name = ['Nanoplatlets','Bulk','Nanocrystalline'] 
#ydata = [nano_pmax_rteo,bulk_pmax_rteo,nanoplt_pmax_rteo]
ydata = [nanoplt_pmax_rteo,bulk_pmax_rteo,nano_pmax_rteo]
xdata = (rteo-rtei)*1e3
color_idx = np.linspace(0.25,0.75,len(mat_name))
for i in range(len(mat_name)):
    plt.plot(xdata, ydata[i]*0.1, color = plt.cm.hot(color_idx[i]), label=mat_name[i], lw=2)
plt.xlabel("TE radial leg length [$mm$]", fontsize=16)
plt.ylabel("Maximum Power Density [$mW/cm^2$]", fontsize=16)
plt.ylim([0,nano_pmax_rteo.max()*0.1*1.05])
plt.xscale('log')
plt.legend(loc=0)	

plt.figure(2)
mat_name1 = ['Nanoplatlets','Bulk','Nanocrystalline'] 
ydata1 = [nanoplt_maxeta_rteo,bulk_maxeta_rteo,nano_maxeta_rteo]
color_idx = np.linspace(0.25,0.75,len(mat_name))
for i in range(len(mat_name)):
    plt.plot(xdata, ydata1[i]*1e2, color = plt.cm.hot(color_idx[i]), label=mat_name[i], lw=2)
plt.xlabel("TE radial leg length [$mm$]", fontsize=16)
plt.ylabel("Maximum Efficiency [%]", fontsize=16)
#plt.ylim([0,nano_eta.max()*1.05])
plt.xscale('log')
plt.legend(loc=0)	

plt.figure(3)
mat_name1 = ['Nanoplatlets','Bulk','Nanocrystalline'] 
ydata2 = [nanoplt_maxp_Rte,bulk_maxp_Rte,nano_maxp_Rte]
ydata3 = [1.0/(nanoplt_hc_maxp*2.0*m.pi*rteo*tw),1.0/(bulk_hc_maxp*2.0*m.pi*rteo*tw),1.0/(nano_hc_maxp*2.0*m.pi*rteo*tw)]
color_idx = np.linspace(0.25,0.75,len(mat_name))
for i in range(len(mat_name)):
    plt.plot(xdata, ydata2[i], color = plt.cm.hot(color_idx[i]), label=mat_name[i], lw=2)
    plt.plot(xdata, ydata3[i],'--', color = plt.cm.hot(color_idx[i]), lw=2)
plt.xlabel("TE radial leg length [$mm$]", fontsize=16)
plt.ylabel("Thermal resistance [$K/W$]", fontsize=16)
#plt.ylim([0,nano_eta.max()*1.05])
plt.xscale('log')
plt.yscale('log')
plt.legend(loc=0)	

plt.figure(4)
mat_name1 = ['Nanoplatlets','Bulk','Nanocrystalline'] 
ydata4 = [nano_maxp_Rte,bulk_maxp_Rte,nanoplt_maxp_Rte]
color_idx = np.linspace(0.25,0.75,len(mat_name))
for i in range(len(mat_name)):
    plt.plot(xdata, ydata4[i], color = plt.cm.hot(color_idx[i]), label=mat_name[i], lw=2)
plt.xlabel("TE radial leg length [$mm$]", fontsize=16)
plt.ylabel("Electrical resistance [$Ohms$]", fontsize=16)
#plt.ylim([0,nano_eta.max()*1.05])
plt.xscale('log')
plt.legend(loc=0)	


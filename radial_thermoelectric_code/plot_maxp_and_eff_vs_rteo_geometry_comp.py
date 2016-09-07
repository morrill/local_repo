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
import seaborn as sns
#%% Material properties:
# temperature
Th = 372.0         # hot stream tempertaure [K]
Tc = 294.0         # cold stream temperature [K]
# geometry
rwi = 3e-3         # inner pipe radius [m]
rwo = 4e-3         # outer pipe
tc = 100e-6                      #contact width
ts = 100e-6                      # spacer thickness for both radial and axial spacers [m]
tw = 10e-3                       # pipe thickness [m]
ti = 100e-6                      # insulator thickness [m]
rtei = rwo + ts                  #inner TE radius [m]
Lw = 2e-3                        # thickness of pipe for radial device and base plate for planar [m]
D = 8e-3                         # outer diemeter of pipe we're copmaring to [m]
FF = np.linspace(0.01,1,num=100) # fill factor array
Ah = tw*m.pi*D                   # hot side area [m^2]
Al = np.zeros(len(FF))
Al = (FF*(Ah/2))              # area of TE legs [m^2]
L = Ah/(2*tw + 2*m.pi*D)       #flat-plate convection characteristic length

# material properties
kw = 385.0                  # pipe thermal conductivity [W/m.K]
ks = 0.1                # spacer/insulator thermal conductivity [W/m.K]
Rc = 1e-9                    # contact resistance [Ohm.m^2]
kins = 0.1                    # insulator thermal conductivity [W/m.K]

# thermal resistances for first two nodes:
Rw = m.log(rwo/rwi)/(2.0*m.pi*tw*kw)   # inner pipe thermal resistance
Rs = m.log(rtei/rwo)/(2.0*m.pi*tw*ks)
Rwpluss = Rw + Rs

Rw = Lw/(Ah*kw)   # inner pipe thermal resistance
Ri = ti/(Ah*kins)
Rwplusi = Rw + Ri
#%% bulk material properties
kn = 1.30                 # n-type material thermal conductivity [W/m.K]
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

# need to run this and raname maxeta_rteo, to start with "material name"_pmax_rteo etc.
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
radial_Qh = (Th-r_T1)/Rwpluss
planar_Qh = (Th-p_T1)/Rwplusi

radial_eta = (radial_p*2.0*m.pi*rtei*tw)/radial_Qh
planar_eta = (planar_p*Ah)/planar_Qh
#%% Get thermal and electrical resistances

Rte = (((2.0*m.pi)/(np.log(rteo/rtei)))*(kn*tn_maxp+ks*ts+kp*tp_maxp))**(-1.0) # TE/spacer total thermal resistance
Relec = ((np.log(rteo/rtei))/(2.0*m.pi*sigma_n*tn_maxp)) + ((np.log(rteo/rtei))/(2.0*m.pi*sigma_p*tp_maxp)) + (Rc/((2*m.pi*((rtei+tc)**2.0 - (rtei)**2.0)) + (2*m.pi*(((rteo)**2.0) - (rteo-tc)**2.0)))) # n- and p-type total electrical resistance plus contact resitance




#%% plot things

plt.figure(1)
mat_name = ['Planar, FF=0.98','Planar, FF=0.20','Radial']
ydata = [pmax_ll,pmax_ll_FFlow,pmax_rteo]
xdata = leglength*1e3
color_idx = np.linspace(0.25,0.75,len(mat_name))
for i in range(len(mat_name)):
    plt.plot(xdata, ydata[i]*0.1, color = plt.cm.hot(color_idx[i]), label=mat_name[i], lw=2)
plt.xlabel("TE leg length [$mm$]", fontsize=16)
plt.ylabel("Maximum Power Density [$mW/cm^2$]", fontsize=16)
plt.ylim([0,pmax_rteo.max()*0.1*1.05])
plt.xscale('log')
plt.legend(loc=0)	

#plt.figure(2)
#mat_name1 = ['Nanoplatlets','Bulk','Nanocrystalline'] 
#ydata1 = [nanoplt_maxeta_rteo,bulk_maxeta_rteo,nano_maxeta_rteo]
#color_idx = np.linspace(0.25,0.75,len(mat_name))
#for i in range(len(mat_name)):
    #plt.plot(xdata, ydata1[i], color = plt.cm.hot(color_idx[i]), label=mat_name[i], lw=2)
#plt.xlabel("TE radial leg length [$mm$]", fontsize=16)
#plt.ylabel("Maximum Efficiency [%]", fontsize=16)
#plt.ylim([0,nano_eta.max()*1.05])
#plt.xscale('log')
#plt.legend(loc=0)	

plt.figure(3)
mat_name1 = ['Planar, FF=0.98','Planar, FF=0.20','Radial'] 
ydata2 = [Rte_FFmax,Rte_FFlow,Rte_Pdmax]
ydata3 = [Rhc_FFmax,Rhc_FFlow,Rhc_Pdmax]
color_idx = np.linspace(0.25,0.75,len(mat_name))
for i in range(len(mat_name)):
    plt.plot(xdata, ydata2[i], color = plt.cm.hot(color_idx[i]), label=mat_name[i], lw=2)
    plt.plot(xdata, ydata3[i],'--', color = plt.cm.hot(color_idx[i]), lw=2)
plt.xlabel("TE leg length [$mm$]", fontsize=16)
plt.ylabel("Thermal resistance [$K/W$]", fontsize=16)
#plt.ylim([0,nano_eta.max()*1.05])
plt.xscale('log')
plt.yscale('log')
plt.legend(loc=0)	

plt.figure(4)
mat_name1 = ['Planar','Planar, FF=0.20','Radial'] 
ydata4 = [Relec_FFmax,Relec_FFlow,Relec_Pdmax]
color_idx = np.linspace(0.25,0.75,len(mat_name))
for i in range(len(mat_name)):
    plt.plot(xdata, ydata4[i]*1e3, color = plt.cm.hot(color_idx[i]), label=mat_name[i], lw=2)
plt.xlabel("TE leg length [$mm$]", fontsize=16)
plt.ylabel("Electrical resistance [$mOhms$]", fontsize=16)
#plt.ylim([0,nano_eta.max()*1.05])
plt.xscale('log')
plt.yscale('log')
plt.legend(loc=0)	


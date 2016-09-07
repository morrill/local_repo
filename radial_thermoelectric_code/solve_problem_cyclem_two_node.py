# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 22:26:26 2016
@author: Michael Orrill
"""
#%% Libraries
from __future__ import division
from scipy.optimize import fsolve
import math as m
import numpy as np
import os.path
import itertools
import csv
import matplotlib.pyplot as plt
import seaborn as sns
#%% 
#water = read_data('water.csv')
air = read_data('air.csv')
# cold side fluid properties (manual input)
#inner_fluid = water           # enter "air" or "water" (w/o "") for the fluid you're working with
outer_fluid = air             # enter "air" or "water" (w/o "") for the fluid you're working with
condition = 'free'       # enter 'forced' or 'free' for forced or free convection, respectively
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
#kh = get_props_for_hh(inner_fluid,Th) # thermal conductivity of fluid in hot side pipe [W/m.K]
kw = 385.0                  # pipe thermal conductivity [W/m.K]
kn = 1.30                  # n-type material thermal conductivity [W/m.K]
kp = 1.35                  # p-type material thermal conductivity [W/m.K]
ks = 0.1                # spacer/insulator thermal conductivity [W/m.K]
Sn = -150e-6                 # n-type Seebeck coefficient [V/K] 
Sp = 225e-6                  # p-type Seebeck coefficient [V/K] 
sigma_n = 86975             # n-type electrical conductivity [S/m]
sigma_p = 90000             # p-type electrical conductivity [S/m]
Rc = 1e-9                    # contact resistance [Ohm.m^2]
# flow conditions
#hh = 4.36*(kh/(2.0*rwi))   # hot side convection coefficient [W/m^2.K] = 4.36*(k/D) 
free_str_vel = 1.0        # free stream velocity [m/s]
# thermal resistances for first two nodes:
Rw = m.log(rwo/rwi)/(2.0*m.pi*tw*kw)   # inner pipe thermal resistance
Rs = m.log(rtei/rwo)/(2.0*m.pi*tw*ks)
Rwpluss = Rw + Rs
#%% solve problem
#rteo = np.linspace(200e-6 + rtei, 1, num=100, dtype=float) # outer TE radius [m] sent in as array of outer radi
elements1 = 32
elements2 = 36
elements3 = 18
elements4 = 90
elements5 = 100
elements6 = 17
elemnetstn = 95

one = np.linspace(200e-6 + rtei, 1e-3 + rtei, elements1, dtype=float, endpoint=False)
two = np.linspace(1e-3 + rtei, 10e-3 + rtei, elements2, dtype=float, endpoint=False)
three = np.linspace(10e-3 + rtei, 100e-3 + rtei, elements3, dtype=float, endpoint=False)
four = np.linspace(100e-3 + rtei, 1000e-3 + rtei, elements4, dtype=float, endpoint=False)
five = np.linspace(1000e-3 + rtei, 2000e-3 + rtei, elements5, dtype=float, endpoint=False)
six = np.linspace(2000e-3 + rtei, 10000e-3 + rtei, elements6, dtype=float, endpoint=True)
rteo = np.concatenate((one,two,three,four,five,six), axis=0)

tn = np.linspace(0.025*tw, 0.975*tw, elemnetstn, dtype=float, endpoint=False) # length of n-type disk along cylinder
mm = np.arange(0.05,2.55,0.05,dtype=float) # ration of electrical load resistance to TE module resistance
tp = (tw-ts)-tn # p-typme material length along cylinder
p = np.zeros((len(mm),len(rteo),len(tn)), dtype=float) # power
rteitemp = np.zeros((len(mm),len(rteo),len(tn)), dtype=float) # outside pipe/inner TE material temperature
rteotemp = np.zeros((len(mm),len(rteo),len(tn)), dtype=float) # outside TE material temperature
outconv = np.zeros((len(mm),len(rteo),len(tn)), dtype=float)  # cold side convection coefficient
pda = np.zeros((len(mm),len(rteo),len(tn)), dtype=float)  # power density(areal) traditional TE definition [W/m^2]
pd = np.zeros((len(mm),len(rteo),len(tn)), dtype=float)  # power density [W/m^3]
eta = np.zeros((len(mm),len(rteo),len(tn)), dtype=float)  # efficiency [%]

# These loop through the problem creating a 3D array for power, unknown temperatures, and outer convection coefficient
# first index (slice) is the mm number, this gives a 2D array of valuse with rows = rteo and columns = tn
for k in range(len(mm)):  
    prog1 = k/len(mm)*100.0
    print 'Progress: %.3f' %prog1
    for j in range(len(tn)):   
        #prog1 = k/len(mm)*100.0+j/len(tn)
        
        for i in range(len(rteo)):
            Rte = (((2.0*m.pi)/(m.log(rteo[i]/rtei)))*(kn*tn[j]+ks*ts+kp*tp[j]))**(-1.0) # TE/spacer total thermal resistance
            Relec = ((m.log(rteo[i]/rtei))/(2.0*m.pi*sigma_n*tn[j])) + ((m.log(rteo[i]/rtei))/(2.0*m.pi*sigma_p*tp[j])) + (Rc/((2*m.pi*((rtei+tc)**2.0 - (rtei)**2.0)) + (2*m.pi*(((rteo[i])**2.0) - (rteo[i]-tc)**2.0)))) # n- and p-type total electrical resistance plus contact resitance
            [hc, T1, T2] = calc_temps_hc(1,rteo[i],Rte,mm[k])
            p[k,i,j] = (mm[k]/(mm[k]+1.0)**2.0)*(Sp-Sn)**2.0*(T1-T2)**2.0/(Relec) # power [W]
            pda[k,i,j] = p[k,i,j]/(2*m.pi*rteo[i]*tw) # power density(areal) traditional TE definition [W/m^2]
            pd[k,i,j] = p[k,i,j]/(m.pi*(rteo[i]**2 - rtei**2)*tw) # power density [W/m^3]
            rteitemp[k,i,j] = T1 # inner TE material temperatures
            rteotemp[k,i,j] = T2 # outer TE material temperatures
            outconv[k,i,j] = hc #outer heat transfer coefficient
            eta[k,i,j] = (p[k,i,j]/((Th - T1)/Rwpluss))*100 # efficiency [%]
            
rteo_maxp_mm = np.zeros(len(mm))
tn_maxp_mm = np.zeros(len(mm))
m_maxp = np.zeros(len(mm))
for i in range(len(mm)):
    [rteo_maxp_mm[i],tn_maxp_mm[i]] = np.unravel_index(p[i].argmax(), p[i].shape) # gets rteo and tn for max p in each mm slice
    m_maxp[i] = p[i].max() # gets max p value from each mm value
maxp = m_maxp.max()        # overall max power
maxpmm = mm[int(m_maxp.argmax())]
maxprteo = rteo[int(rteo_maxp_mm[m_maxp.argmax()])]
maxptn = tn[int(tn_maxp_mm[m_maxp.argmax()])]

T1atmax = rteitemp[int(m_maxp.argmax()),int(rteo_maxp_mm[m_maxp.argmax()]),int(tn_maxp_mm[m_maxp.argmax()])]
T2atmax = rteotemp[int(m_maxp.argmax()),int(rteo_maxp_mm[m_maxp.argmax()]),int(tn_maxp_mm[m_maxp.argmax()])]
Tm = np.mean([T1atmax,T2atmax])
print 'Max power density = %.3f [W/m^2]\nT1 = %.3f [K]\nT2 = %.3f [K]\nAverage temp = %.3f [K]\ntn at maxp = %.3f [mm]\nrteo at maxp = %.3f [mm]\nmm at maxp = %.3f' %(p.max(),T1atmax,T2atmax,Tm,maxptn*1e3,maxprteo*1e3,maxpmm)

plt.figure(1)
plt.xlabel('m',fontsize=16)
plt.ylabel('Power density [$W/m^2$]',fontsize=16)
plt.plot(mm,m_maxp,'r',lw=2)
#%% plots for sliced data to give more accurate curves of of power density vs. leg length by finding
#mm values and tn values that maximizes power density for each value of rteo
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

print 'Max power density = %.3f [W/m^2]\nT2 = %.3f [K]\nAverage temp = %.3f [K]\ntn at maxp = %.3f [mm]\nrteo at maxp = %.3f [mm]\nmm at maxp = %.3f' %(p.max(),T2atmax,Tm,maxptn2*1e3,maxprteo2*1e3,maxpmm2)

#plot contour of max power density vs outer radius and n-type material thickness and a line graph of max 
#power density vs rteo
tnplot = tn/tw
[X,Y] = np.meshgrid(tnplot,(rteo-rtei)*1e3)
#size = 10
plt.figure(2)
#plt.subplot(121)
plt.ylim((rteo[0] - rtei)*1e3, (rteo[-1] - rtei)*1e3)
plt.xlim(tnplot[0], tnplot[-1])
plt.xlabel('n-type material thickness/total device thickness',fontsize=16)
plt.ylabel('Leg Length [$mm$]',fontsize=16)
# flips axis because of how p is defined with rows (y axis) as rteo and columns (x axis) as tn
contf = plt.contourf(X,Y, pmax_rteo_tn*0.1, levels=np.linspace(0, int(round(maxp2*0.1)), 25), extend='both', cmap = sns.light_palette((18, 94, 50), input="husl", as_cmap=True))
cbar = plt.colorbar(contf)
plt.title('Power Density [$mW/cm^2$]', fontsize=16)
cbar.set_label('Power Density [$mW/cm^2$]', fontsize=16)
cbar.set_ticks(np.linspace(0, int(maxp*0.1), 10));
plt.yscale('log')

plt.figure(3)
plt.xlabel('Leg Length [$mm$]',fontsize=16)
plt.ylabel('Max Power Density [$mW/cm^2$]',fontsize=16)
plt.plot((rteo - rtei)*1e3, pmax_rteo*0.1,'r',lw=2)
plt.xscale('log')

plt.show()
#%% more output data
#thermal reistance values at max power density
Rte = (((2.0*m.pi)/(m.log(rteo[int(rteo_maxp_mm[m_maxp.argmax()])]/rtei)))*(kn*tn[int(tn_maxp_mm[m_maxp.argmax()])]+ks*ts+kp*tp[int(tn_maxp_mm[m_maxp.argmax()])]))**(-1.0)
Rw = m.log(rtei/rwi)/(2.0*m.pi*tw*kw)
#Rhh = (hh*2.0*m.pi*rwi*tw)**-1.0
hc = outconv[int(m_maxp.argmax()),int(rteo_maxp_mm[m_maxp.argmax()]),int(tn_maxp_mm[m_maxp.argmax()])]
Rhc = (hc*2.0*m.pi*rteo[int(rteo_maxp_mm[m_maxp.argmax()])]*tw)**-1.0
print 'Rte = %.3f [W/K]\nRw = %.3f [W/K]\nhc = %.3f [W/m^2K]\nRhc = %.3f [W/K]' %(Rte,Rw,hc,Rhc)

#electrical values
Relec_te = ((m.log(rteo[int(rteo_maxp_mm[m_maxp.argmax()])]/rtei))/(2.0*m.pi*sigma_n*tn[int(tn_maxp_mm[m_maxp.argmax()])])) + ((m.log(rteo[int(rteo_maxp_mm[m_maxp.argmax()])]/rtei))/(2.0*m.pi*sigma_p*tp[int(tn_maxp_mm[m_maxp.argmax()])]))
Relec_contact = (Rc/((2*m.pi*((rtei+tc)**2.0 - (rtei)**2.0)) + (2*m.pi*(((rteo[int(rteo_maxp_mm[m_maxp.argmax()])])**2.0) - (rteo[int(rteo_maxp_mm[m_maxp.argmax()])]-tc)**2.0))))
print 'Relec_te = %.3f [mOhms]\nRelec_contact = %.3f [microOhms]\nRelec = %.3f [mOhms]' %(Relec_te*1e3, Relec_contact*1e6,Relec*1e3)

#module ZTm
ZTm = (Sp - Sn)**2*Tm/(Relec*Rte)
#%% Heat flow model check
#Check model
#Qh = (Th - Twatmax)/Rhh
Qh = (Th - T1atmax)/Rwpluss
Qte = (T1atmax - T2atmax)/Rte
Half_I2Re = (Sp-Sn)**2.0*(T1atmax-T2atmax)**2.0/(2.0*Relec*(mm[int(m_maxp.argmax())]+1.0)**2.0)
SIT1 = (Sp-Sn)**2.0*(T1atmax-T2atmax)*T1atmax/(Relec*(mm[int(m_maxp.argmax())]+1.0))
SIT2 = (Sp-Sn)**2.0*(T1atmax-T2atmax)*T2atmax/(Relec*(mm[int(m_maxp.argmax())]+1.0))
Qte2 = Qte - Half_I2Re + SIT1
Qc = (T2atmax - Tc)/Rhc
P = Qh - Qc
Pd = P/(2.0*m.pi*rtei*tw)

print 'Qh = %.6f [W]\nQte - 1/2I^2Re + SIT1 = %.6f [W]\nQc = %.6f [W]\nP = %.6f [W]\nPd = %.6f [W/m^2]' %(Qh,Qte2,Qc,P, Pd)
#%% All user functions
def calc_temps_hc(hc,rteo,Rte,mm):
    D = 2.0*rteo
    # average temperature for initial guessing purposes
    Tguess = (np.mean([Th,Tc]),np.mean([Th,Tc]))  
    #count = 0     
    while True:   
        [T1, T2] = solve_temps(hc,rteo,Rte,mm,Tguess)
        Tguess = (T1,T2)
        #print "Tw = %0.3f\nT1 = %0.3f\nT2 = %0.3f" %(Tw,T1,T2)
        hc_new = h_bar(T2,D) 
        t = abs(hc - hc_new)
        #print t
        #count += 1.0
        if t < 1e-6:
            #print 'iterations: %0f' %count
            break
        hc = hc_new
    #print "Tw = %0.3f\nT1 = %0.3f\nT2 = %0.3f" %(Tw,T1,T2)    
    #print "hc = %0.3f" %hc_new
    return hc_new,T1,T2
 
def solve_temps(hc,rteo,Rte,mm,Tguess):
# temperature solver
    def f(p):
    # funtion for equations    
        T1,T2 = p
        return ((1/Rwpluss)*(Th-T1) - (T1-T2)/Rte + (Sp-Sn)**2.0*(T1-T2)**2.0/(2.0*Relec*(mm+1.0)**2.0) - (Sp-Sn)**2.0*(T1-T2)*T1/(Relec*(mm+1.0)), (T1-T2)/Rte + (Sp-Sn)**2.0*(T1-T2)**2.0/(2.0*Relec*(mm+1.0)**2.0) + (Sp-Sn)**2.0*(T1-T2)*T2/(Relec*(mm+1.0)) - hc*2.0*m.pi*rteo*tw*(T2-Tc))    
    T1,T2 = fsolve(f,Tguess)
    #print "Method 1:\nTw = %0.3f\nT1 = %0.3f\nT2 = %0.3f" %(Tw,T1,T2)
    return T1,T2
    
def h_bar(T2,D):
# gets h_bar for iterative problem
    T_f = np.mean([Tc,T2]) 
    [visc, k, Prs, Prf, Pr] = get_props(outer_fluid, Tc, T2, T_f)
    Re = free_str_vel*D/visc
    if condition == 'forced':
        hC = Churchill_Burnstein(Prf, Re)*(k/D)
        return hC
    elif condition == 'free':
        hC = Churchill_Chu(T_f, T2, D, visc, Prf)*(k/D)
        return hC

def get_props(medium, Tc, T2, T_f):
# fluid properties from thermophysical tables
    data = medium
    T = data.keys()
    T_val = []
    for i in range(len(T)):
        T_val.append(int(T[i]))
    T_val.sort()
    if not min(T_val) <= T2 <= max(T_val):                    # check if surface temperature is in tables range
        print "Error: Surface temperature, T2 = %s is not in range of tabular data, choose temperature in the range of %s - %s" %(Tc,min(T_val),max(T_val))
        return
    else:
        if str(T2) in T and str(Tc) in T and str(T_f) in T: # no interpolation needed
            visc = float(data[str(int(T2))]['visc'])         # kinimatic viscosity 
            k = float(data[str(int(T2))]['thermcon'])        # thermal conductivity
            Prs = float(data[str(int(T2))]['Pr'])            # Prandtl number at surface
            Pr = float(data[str(int(Tc))]['Pr'])             # Prandtl number in free stream
            Prf = float(data[str(int(T_f))]['Pr'])           # Prandtl number at film temp
            return visc, k, Prs, Prf, Pr
        else: # interpolation needed
            for lower,upper in zip (T_val[:-1],T_val[1:]): # get value above and below T2
                if lower <= T2 <= upper:
                    [l2,u2] = lower, upper
                    break
            for lower,upper in zip (T_val[:-1],T_val[1:]): # get value above and below T_f
                if lower <= T_f <= upper:
                    [lf,uf] = lower, upper
                    break
            for lower,upper in zip (T_val[:-1],T_val[1:]): # get value above and below Tc
                if lower <= Tc <= upper:
                    [lc,uc] = lower, upper
                    break
            # assign fluid properties
            visc = float(data[str(lc)]['visc']) + (float(data[str(uc)]['visc'])-float(data[str(lc)]['visc']))*((Tc-lc)/(uc-lc))
            k = float(data[str(lc)]['thermcon']) + (float(data[str(uc)]['thermcon'])-float(data[str(lc)]['thermcon']))*((Tc-lc)/(uc-lc))
            Pr = float(data[str(lc)]['Pr']) + (float(data[str(uc)]['Pr'])-float(data[str(lc)]['Pr']))*((Tc-lc)/(uc-lc))
            Prs = float(data[str(l2)]['Pr']) + (float(data[str(u2)]['Pr'])-float(data[str(l2)]['Pr']))*((T2-l2)/(u2-l2))
            Prf = float(data[str(lf)]['Pr']) + (float(data[str(uf)]['Pr'])-float(data[str(lf)]['Pr']))*((T_f-lf)/(uf-lf))
            return visc, k, Prs, Prf, Pr
         
def get_props_for_hh(medium, Th):
# fluid properties from thermophysical tables
    data = medium
    T = data.keys()
    T_val = []
    for i in range(len(T)):
        T_val.append(int(T[i]))
    T_val.sort()
    if str(Th) in T:                          # no interpolation needed
        k = float(data[str(int(Th))]['thermcon'])  # thermal conductivity
        return k
    else: # interpolation needed
        for lower,upper in zip (T_val[:-1],T_val[1:]): # get value above and below T_s
            if lower <= Th <= upper:
                [l,u] = lower, upper
                break
        # assign fluid properties
        k = float(data[str(l)]['thermcon']) + (float(data[str(u)]['thermcon'])-float(data[str(l)]['thermcon']))*((Th-l)/(u-l))
        return k
          
def read_data(filename):
# Reads in data for air or water
# outputs dictionary with temperatures as keys that reference viscosity, thermal conductivity, and Prandtl #
    path = os.path.join(os.getcwd(),filename) 
    if os.path.isfile(path):
        infile = open(path)
        row_num = len(infile.readlines())
        infile = open(path)
        f = csv.reader(infile)
        b = list(f)
        tempK = []
        for i in range(1,len(b)):
            tempK.append(b[i][0])
        count = 0
        datahead = ['visc','thermcon','Pr']
        prop = {head:{i:[] for i in datahead} for head in tempK}
        infile = open(path)
        f = csv.reader(infile)
        for row in itertools.islice(f,1,row_num):
            prop[tempK[count]]['visc'] = row[1]
            prop[tempK[count]]['thermcon'] = row[2]
            prop[tempK[count]]['Pr'] = row[3]
            count += 1
        infile.close()
        return prop
    else:
        print "File: "+ file[59:] +" not found in directory"
    return
    
def Churchill_Burnstein(Prf, Re):
# Churchill and Burnstein method
# all properties at film temperature
    Nu = 0.3+(0.62*Re**(1/2)*Prf**(1/3))/(1+(0.4/Prf)**(2/3))**(1/4)*(1+(Re/282000)**(5/8))**(4/5)
    return Nu

# Free convection correlations
def Churchill_Chu(T_f, T2, D, visc, Prf):
# Churchill and Chu method
    Ra = 9.81*(1/T_f)*(T2-Tc)*D**3/(visc**2)*Prf
    Nu = (0.60 + 0.387*Ra**(1/6)/(1+(0.559/Prf)**(9/16))**(8/27))**2
    return Nu
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 22:26:26 2016
The next thing to do is to make rtei into an array of values and solve temperature and conv. coeff. for all 
so that power may be computed.
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
#%% 
# cold side fluid properties (manual input)
inner_fluid = 'water'             # enter 'air' or 'water' for the fluid you're working with
outer_fluid = 'air'             # enter 'air' or 'water' for the fluid you're working with
condition = 'forced'       # enter 'forced' or 'free' for forced or free convection, respectively
# temperature
Th = 474.0         # hot stream tempertaure [K]
Tc = 274.0         # cold stream temperature [K]
# average temperature for initial guessing purposes
Tguess = (np.mean([Th,Tc]),np.mean([Th,Tc]),np.mean([Th,Tc]))
# geometry
rwi = 0.05          # inner pipe radius [m]
rtei = rwi +0.005   # outer pipe/inner TE radius [m]
tn = 0.4            # n-type thickness [m]
tp = 0.4            # p-type thickness [m]
ts = abs(1.0-tn-tp) # spacer thickness [m]
tw = tn+ts+tp       # pipe thickness [m]
#D = 2*rteo          # outer diameter of TE tube
# material properties
kh = get_props_for_hh(inner_fluid,Th) # thermal conductivity of fluid in hot side pipe [W/m.K]
kw = 1.0                  # pipe thermal conductivity [W/m.K]
kn = 1.0                  # n-type material thermal conductivity [W/m.K]
kp = 1.0                  # p-type material thermal conductivity [W/m.K]
ks = 1.0                  # spacer/insulator thermal conductivity [W/m.K]
Sn = -1.0                 # n-type Seebeck coefficient [µV/K] 
Sp = 1.0                  # p-type Seebeck coefficient [µV/K] 
sigma_n = 1.0             # n-type electrical conductivity []
sigma_p = 1.0             # p-type electrical conductivity []
# flow conditions
hh = 4.36*(kh/(2.0*rwi))   # hot side convection coefficient [W/m^2.K] = 4.36*(k/D) 
free_str_vel = 1.0        # free stream velocity [m/s]
# thermal resistances:
#Rte = (m.log(rteo/rtei)/(2.0*m.pi))*(1.0/(kn*tn)+1.0/(ks*ts)+1.0/(kp*tp)) # TE/spacer total R
Rw = m.log(rtei/rwi)/(2.0*m.pi*tw*kw)                                    # inner pipe thermal resistance
# electrical resistance
#Relec = (m.log(rteo/rtei)/(2.0*m.pi))*(1.0/(sigma_n*tn)+1.0/(sigma_p*tp))    # n- and p-type total Re
#%% solve problem
rteo = np.arange(1.1*rtei, 1e2, rtei, dtype=float) # outer TE radius [m] sent in as array of outer radi MAKE AN ARRAY EVENTUALLY
p = np.zeros(len(rteo))
walltemp = np.zeros(len(rteo))
rteitemp = np.zeros(len(rteo))
rteotemp = np.zeros(len(rteo))
deltaT = np.zeros(len(rteo))
outconv = np.zeros(len(rteo))
for i in range(len(rteo)):
    Rte = (m.log(rteo[i]/rtei)/(2.0*m.pi))*(1.0/(kn*tn)+1.0/(ks*ts)+1.0/(kp*tp)) # TE/spacer total thermal resistance
    Relec = (m.log(rteo[i]/rtei)/(2.0*m.pi))*(1.0/(sigma_n*tn)+1.0/(sigma_p*tp)) # n- and p-type total electrical resistance
    [hc, Tw, T1, T2] = calc_temps_hc(1,rteo[i],Rte)
    p[i] = (Sp-Sn)**2.0*(T1-T2)**2.0/Relec
    walltemp[i] = Tw
    rteitemp[i] = T1
    rteotemp[i] = T2
    deltaT[i] = abs(T1-T2)
    outconv[i] = hc
#%% plot results
plt.plot(rteo,p) # power vs rteo
plt.title("Power vs. outer radius")
plt.xlabel("Outer radius [m]")
plt.ylabel("Power [W]")

plt.figure()
plt.plot(rteo,walltemp, label='Pipe wall') # inner pipe wall temp vs rteo
plt.xlabel("Outer radius [m]")
plt.ylabel("Inner pipe temperature [K]")

plt.plot(rteo,rteitemp, label='T1') # hot side TE temp vs rteo
plt.xlabel("Outer radius [m]")
plt.ylabel("Inner (hot side) TE temperature [K]")

plt.plot(rteo,rteotemp, label='T2') # cold side TE temp vs rteo
plt.xlabel("Outer radius [m]")
plt.ylabel("Outer (cold side) TE temperature [K]")

plt.plot(rteo,deltaT, label='T1-T2') # cold side TE temp vs rteo
plt.xlabel("Outer radius [m]")
plt.ylabel("TE temperature difference [K]")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0.)

plt.figure()
plt.title("Hot side convection coefficient vs. outer radius")
plt.plot(rteo,outconv) # outer convection coefficient vs rteo
plt.xlabel("Outer radius [m]")
plt.ylabel("Cold side convection coefficient [W/m^2.K]")
#%% All user functions
def calc_temps_hc(hc,rteo,Rte):
    D = 2.0*rteo       
    while True:   
        [Tw, T1, T2] = solve_temps(hc,rteo,Rte)
        #print "Tw = %0.3f\nT1 = %0.3f\nT2 = %0.3f" %(Tw,T1,T2)
        hc_new = h_bar(T2,D) 
        t = abs(hc - hc_new)
        #print t
        if t < 1e-5:
            break
        hc = hc_new
    #print "Tw = %0.3f\nT1 = %0.3f\nT2 = %0.3f" %(Tw,T1,T2)    
    #print "hc = %0.3f" %hc_new
    return hc_new, Tw, T1, T2
 
def solve_temps(hc,rteo,Rte):
# temperature solver
    def f(p):
    # funtion for equations    
        Tw,T1,T2 = p
        return ((hh*2.0*m.pi*rwi*tw)*(Th-Tw) - 2.0*m.pi*kw*tw/m.log(rtei/rwi)*(Tw-T1), 2.0*m.pi*kw*tw/m.log(rtei/rwi)*(Tw-T1) - (T1-T2)/Rte + (Sp-Sn)**2.0*(T1-T2)**2.0/(4.0*Relec) - (Sp-Sn)**2.0*(T1-T2)*T1/(2.0*Relec), (T1-T2)/Rte + (Sp-Sn)**2.0*(T1-T2)**2.0/(4.0*Relec) + (Sp-Sn)**2.0*(T1-T2)*T2/(2.0*Relec) - hc*2.0*m.pi*rteo*tw*(T2-Tc))    
    Tw,T1,T2 = fsolve(f,Tguess)
    #print "Method 1:\nTw = %0.3f\nT1 = %0.3f\nT2 = %0.3f" %(Tw,T1,T2)
    return Tw, T1, T2
    
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
    data = read_data(medium + '.csv')
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
    data = read_data(medium + '.csv')
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
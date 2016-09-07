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
#%% 
# read in fluid property data:
water = read_data('water.csv')
air = read_data('air.csv')
# cold side fluid properties (manual input)
inner_fluid = water           # enter "air" or "water" (w/o "") for the fluid you're working with
outer_fluid = air             # enter "air" or "water" (w/o "") for the fluid you're working with
condition = 'free'       # enter 'forced' or 'free' for forced or free convection, respectively
# temperature
Th = 372.0         # hot stream tempertaure [K]
Tc = 292.0         # cold stream temperature [K]
# geometry
rwi = 12e-3/2       # inner pipe radius [m]
rtei = 14e-3/2       # outer pipe/inner TE radius [m]
#tn = 40e-6          # n-type thickness [m]
#tp = 40e-6          # p-type thickness [m]
ts = 25e-6          # spacer thickness [m]
tw = 10000e-6 #tn+ts+tp       # pipe thickness [m]
#D = 2*rteo         # outer diameter of TE tube
# material properties
mm = 1.13                     # ratio of load resistance to TE resistance
kh = get_props_for_hh(inner_fluid,Th) # thermal conductivity of fluid in hot side pipe [W/m.K]
kw = 385.0                  # pipe thermal conductivity [W/m.K]
kn = 1.30                  # n-type material thermal conductivity [W/m.K]
kp = 1.35                  # p-type material thermal conductivity [W/m.K]
ks = 1.0e-1                  # spacer/insulator thermal conductivity [W/m.K]
Sn = -150e-6                 # n-type Seebeck coefficient [V/K] 
Sp = 225e-6                  # p-type Seebeck coefficient [V/K] 
sigma_n = 86957             # n-type electrical conductivity [S/m]
sigma_p = 90000             # p-type electrical conductivity [S/m]
# flow conditions
hh = 4.36*(kh/(2.0*rwi))   # hot side convection coefficient [W/m^2.K] = 4.36*(k/D) 
free_str_vel = 1.0        # free stream velocity [m/s]
# thermal resistances:
Rw = m.log(rtei/rwi)/(2.0*m.pi*tw*kw)                                    # inner pipe thermal resistance
#%% solve problem
#==============================================================================
# elements = 50
# one = np.linspace(200e-6 + rtei, 1e-3 + rtei, elements, dtype=float)
# two = np.linspace(1e-3 + rtei, 10e-3 + rtei, elements, dtype=float)
# three = np.linspace(10e-3 + rtei, 100e-3 + rtei, elements, dtype=float)
# four = np.linspace(100e-3 + rtei, 1000e-3 + rtei, elements, dtype=float)
# five = np.linspace(1000e-3 + rtei, 10000e-3 + rtei, elements, dtype=float, endpoint=True)
# rteo = np.concatenate((one,two,three,four,five), axis=0)
#==============================================================================
elements1 = 16
elements2 = 18
elements3 = 18
elements4 = 90
elements5 = 100
elements6 = 17
elemnetstn = 100
one = np.linspace(200e-6 + rtei, 1e-3 + rtei, elements1, dtype=float, endpoint=False)
two = np.linspace(1e-3 + rtei, 10e-3 + rtei, elements2, dtype=float, endpoint=False)
three = np.linspace(10e-3 + rtei, 100e-3 + rtei, elements3, dtype=float, endpoint=False)
four = np.linspace(100e-3 + rtei, 1000e-3 + rtei, elements4, dtype=float, endpoint=False)
five = np.linspace(1000e-3 + rtei, 2000e-3 + rtei, elements5, dtype=float, endpoint=False)
six = np.linspace(2000e-3 + rtei, 10000e-3 + rtei, elements6, dtype=float, endpoint=True)
rteo = np.concatenate((one,two,three,four,five,six), axis=0)

tn = np.linspace(0.01*(tw-ts), 0.99*(tw-ts), elemnetstn, dtype=float, endpoint=True) # length of n-type disk along cylinder
mm = np.arange(0.01,2.01,0.01,dtype=float) # ration of electrical load resistance to TE module resistance

tn = np.linspace(0.01*(tw-ts), 0.99*(tw-ts), num=len(rteo), dtype=float)
tp = (tw-ts)-tn
p = np.zeros((len(rteo),len(rteo)), dtype=float)
walltemp = np.zeros((len(rteo),len(rteo)), dtype=float)
rteitemp = np.zeros((len(rteo),len(rteo)), dtype=float)
rteotemp = np.zeros((len(rteo),len(rteo)), dtype=float)
deltaT = np.zeros((len(rteo),len(rteo)), dtype=float)
outconv = np.zeros((len(rteo),len(rteo)), dtype=float)
for j in range(len(tn)):
    prog = j/len(tn)*100.0    
    print 'Progress: %.3f %%' %prog
    for i in range(len(rteo)):
        Rte = ((2.0*m.pi)/(m.log(rteo[i]/rtei))*(kn*tn[j]+ks*ts+kp*tp[j]))**-1.0 # TE/spacer total thermal resistance
        #Relec = (m.log(rteo[i]/rtei)/(2.0*m.pi))*(1.0/(sigma_n*tn[j])+1.0/(sigma_p*tp[j])) # n- and p-type total electrical resistance
        Relec = (m.log(rteo[i]/rtei)/(2.0*m.pi))*(1.0/(sigma_n*tn[j])+1.0/(sigma_p*tp[j])) + ts*(1/(sigma_c*m.pi*((rtei+100e-6)**2.0 - (rtei)**2.0)) + 1/(sigma_c*m.pi*((rteo[i])**2.0 - (rteo[i]-100e-6)**2.0)))
        [hc, Tw, T1, T2] = calc_temps_hc(1,rteo[i],Rte)
        p[i,j] = (mm/(mm+1.0)**2.0)*(Sp-Sn)**2.0*(T1-T2)**2.0/(Relec*2.0*m.pi*rtei*tw) # power density [W/m^2]
        walltemp[i,j] = Tw # pipe inner wall temperatures
        rteitemp[i,j] = T1 # inner TE material temperatures
        rteotemp[i,j] = T2 # outer TE material temperatures
        deltaT[i,j] = abs(T1-T2)
        outconv[i,j] = hc
#%% find max power value and index to get the tn/rteo that generated it
[rteo_maxp,tn_maxp] = np.unravel_index(p.argmax(), p.shape) # assign indicies for rteo and tn at max power
Twatmax = walltemp[rteo_maxp,tn_maxp]
T1atmax = rteitemp[rteo_maxp,tn_maxp]
T2atmax = rteotemp[rteo_maxp,tn_maxp]
Tm = np.mean([Twatmax,T1atmax,T2atmax])
print 'Max power density = %.3f [W/m^2]\nWall temp = %.3f\nT1 = %.3f\nT2 = %.3f\nAverage temp = %.3f' %(p.max(),Twatmax,T1atmax,T2atmax,Tm)

#plot contour of power density vs outer radius and n-type material thickness
tnplot = tn/tw
[X,Y] = np.meshgrid(tnplot,(rteo-rtei)*1e3)
size = 10
plt.figure(1,figsize=(2*size,size))
plt.subplot(121)
plt.ylim((rteo[0]-rtei)*1e3, (rteo[-1]-rtei)*1e3)
plt.xlim(tnplot[0], tnplot[-1])
plt.xlabel('n-type material thickness/total device thickness',fontsize=16)
plt.ylabel('TE material radial thickness [$mm$]',fontsize=16)
# flips axis because of how p is defined with rows (y axis) as rteo and columns (x axis) as tn
contf = plt.contourf(X,Y, p, levels=np.linspace(0, p.max(), 100))
cbar = plt.colorbar(contf)
plt.title(r'Power Density [$W/m^2$]', fontsize=16)
cbar.set_label('Power Density [$W/m^2$]', fontsize=16)
cbar.set_ticks(np.linspace(0, p.max(), 15));
plt.yscale('log')
# max power vs leg length at thickness that maximized power
plt.subplot(122)
plt.xlabel('TE material radial thickness [$mm$]',fontsize=16)
plt.ylabel('Power Density [$W/m^2$]',fontsize=16)
plt.plot((rteo-rtei)*1e3,p[:,tn_maxp],'r')
plt.xlim([rteo.min()*1e3,rteo.max()*1e3])
plt.xscale('log')

#%% plot results
u = np.arange(0,len(rteo),20)
plt.figure(2,figsize=(3*size,size))
plt.subplot(131)
plt.plot(rteo*1e3,p[:,u]*(0.1)) # power vs rteo [mW/cm^2]
plt.xscale('log')
plt.xlim([rteo[0]*1e3,rteo[-1]*1e3])
plt.title("Power Density vs. outer radius")
plt.xlabel("Outer radius [mm]")
plt.ylabel("Power Density [mW/cm^2]")

plt.subplot(132)
plt.xscale('log')
plt.xlabel("Outer radius [m]")
plt.ylabel("Temperature [K]")
plt.plot(rteo,walltemp[:,u], label='Pipe wall') # inner pipe wall temp vs rteo
plt.plot(rteo,rteitemp[:,u], label='T1') # hot side TE temp vs rteo
plt.plot(rteo,rteotemp[:,u], label='T2') # cold side TE temp vs rteo
plt.ylim(Tc,Th)
#plt.plot(rteo,deltaT, label='T1-T2') # cold side TE temp vs rteo
#plt.legend(loc=0)

plt.subplot(133)
plt.xscale('log')
plt.title("Hot side convection coefficient vs. outer radius")
plt.plot(tn,outconv[:,u]) # outer convection coefficient vs rteo
plt.xlabel("Outer radius [m]")
plt.ylabel("Cold side convection coefficient [W/m^2.K]")
#%% All user functions
def calc_temps_hc(hc,rteo,Rte):
    D = 2.0*rteo
    # average temperature for initial guessing purposes
    Tguess = (np.mean([Th,Tc]),np.mean([Th,Tc]),np.mean([Th,Tc]))
    count = 0    
    while True:   
        [Tw, T1, T2] = solve_temps(hc,rteo,Rte,Tguess)
        Tguess = (Tw,T1,T2)
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
    return hc_new, Tw, T1, T2
 
def solve_temps(hc,rteo,Rte,Tguess):
# temperature solver
    def f(p):
    # funtion for equations    
        Tw,T1,T2 = p
        return ((hh*2.0*m.pi*rwi*tw)*(Th-Tw) - 2.0*m.pi*kw*tw/m.log(rtei/rwi)*(Tw-T1), 2.0*m.pi*kw*tw/m.log(rtei/rwi)*(Tw-T1) - (T1-T2)/Rte + (Sp-Sn)**2.0*(T1-T2)**2.0/(2.0*Relec*(mm+1.0)**2.0) - (Sp-Sn)**2.0*(T1-T2)*T1/(Relec*(mm+1.0)), (T1-T2)/Rte + (Sp-Sn)**2.0*(T1-T2)**2.0/(2.0*Relec*(mm+1.0)**2.0) + (Sp-Sn)**2.0*(T1-T2)*T2/(Relec*(mm+1.0)) - hc*2.0*m.pi*rteo*tw*(T2-Tc))    
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
    data = medium #read_data(medium + '.csv')
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
    data = medium #read_data(medium + '.csv')
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
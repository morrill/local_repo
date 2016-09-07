# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 15:31:24 2016
Calculating heat transfer coefficient for radial thermoelectric project.
author: Michael Orrill

Guess T2, plug into equations, solve Tw, T1, and hc, then compare h to the h_bar() calc. here.
"""
#%% import libraries
from __future__ import division
import numpy as np
import os.path
import itertools
import csv
#%% Manually input flow conditions and fluid properties
# fluid properties
visc = 20.92e-6            # kinimatic viscosity 
free_str_vel = 10          # free stream velocity
k = 30e-3                  # thermal conductivity
T_inf = 300                # free stream temperature
T_s = 401                  # surface temperture
T_f = np.mean([T_inf,T_s]) # film temperature
Pr = 0.707                 # Prandtl number @ T_inf
Prs = 0.690                # Prandtl number @ T_s
Prf = 0.7                  # Prandtl number @ T_f
# Geometry of cylinder
rout = 0.1          # outer radius of TE material
D = 12.7e-3         # outer diameter of TE tube

# Reynolds number for flow conditions
Re = free_str_vel*D/visc
print"Reynolds number = %0.2f" %Re
#%% Automatically call for fluid properties given surface and surrounding tempertaure
# fluid properties (manual input)
medium = 'air'             # enter 'air' or 'water' for the fluid you're working with
condition = 'forced'       # enter 'forced' or 'free' for forced or free convection, respectively
T_inf = 100                # free stream temperature [K]
T_s = 146.909                # surface temperture [K]
T_f = np.mean([T_inf,T_s]) # film temperature [K]
free_str_vel = 5           # free stream velocity [m/s]
# Geometry of cylinder
rout = 0.01                # outer radius of TE material [m]
D = 2*rout                 # outer diameter of TE tube

# this is automatic
[visc, k, Prs, Prf, Pr] = get_props(medium, T_inf, T_s, T_f)
Re = free_str_vel*D/visc
h_bar()

#print"Reynolds number = %0.2f\nDynamic viscosity = %0.10f\nThermal conductivity = %0.10f\nPrandtl number in free stream = %0.2f\nPrandtl number at surface = %0.2f\nPrandtl number at film temp = %0.2f\n" %(Re,visc,k,Pr,Prs,Prf)
#%% All user functions
def h_bar():
# gets h_bar for iterative problem
    if condition == 'forced':
        if Re < 0.4 or 4e5 < Re <= 1e6:
            hZ = Zukauskas()*(k/D)
            hC = Churchill_Burnstein()*(k/D)
            return hZ, hC
        elif Re < 1 or Re > 1e6:
            hC = Churchill_Burnstein()*(k/D)
            return hC
        elif 0.4 <= Re <= 4e5:
            hZ = Zukauskas()*(k/D)
            hC = Churchill_Burnstein()*(k/D)
            hH = Hilpert()*(k/D)
            return hZ, hC, hH
    elif condition == 'free':
        hM = Morgan()*(k/D)
        hC = Churchill_Chu()*(k/D)
        return hM, hC

def get_props(medium, T_inf, T_s, T_f):
# fluid properties from thermophysical tables
    data = read_data(medium + '.csv')
    T = data.keys()
    T_val = []
    for i in range(len(T)):
        T_val.append(int(T[i]))
    T_val.sort()
    if not min(T_val) <= T_s <= max(T_val):                    # check if surface temperature is in tables range
        print "Error: Surface temperature, T_s = %s is not in range of tabular data, choose temperature in the range of %s - %s" %(T_s,min(T_val),max(T_val))
        return
    else:
        if str(T_s) in T:                          # no interpolation needed
            visc = float(data[str(T_s)]['visc'])   # kinimatic viscosity 
            k = float(data[str(T_s)]['thermcon'])  # thermal conductivity
            Prs = float(data[str(T_s)]['Pr'])      # Prandtl number at surface
            Pr = float(data[str(T_inf)]['Pr'])     # Prandtl number in free stream
            Prf = float(data[str(T_f)]['Pr'])      # Prandtl number at film temp
            return visc, k, Prs, Prf, Pr
        else: # interpolation needed
            for lower,upper in zip (T_val[:-1],T_val[1:]): # get value above and below T_s
                if lower <= T_s <= upper:
                    [l,u] = lower, upper
                    break
            for lower,upper in zip (T_val[:-1],T_val[1:]): # get value above and below T_f
                if lower <= T_f <= upper:
                    [lf,uf] = lower, upper
                    break
            # assign fluid properties
            visc = float(data[str(l)]['visc']) + (float(data[str(u)]['visc'])-float(data[str(l)]['visc']))*((T_s-l)/(u-l))
            k = float(data[str(l)]['thermcon']) + (float(data[str(u)]['thermcon'])-float(data[str(l)]['thermcon']))*((T_s-l)/(u-l))
            Prs = float(data[str(l)]['Pr']) + (float(data[str(u)]['Pr'])-float(data[str(l)]['Pr']))*((T_s-l)/(u-l))
            Prf = float(data[str(lf)]['Pr']) + (float(data[str(uf)]['Pr'])-float(data[str(lf)]['Pr']))*((T_s-lf)/(uf-lf))
            Pr = float(data[str(T_inf)]['Pr'])
            return visc, k, Prs, Prf, Pr

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
    
# Forced convection correlations
def Hilpert():
# Hilper method
    if 0.4 <= Re < 4:
        C = 0.989
        m = 0.330
    elif 4 <= Re < 40:
        C = 0.911
        m = 0.385
    elif 40 <= Re < 4e3:
        C = 0.683
        m = 0.466
    elif 4e3 <= Re < 4e4:
        C = 0.193
        m = 0.618
    elif 4e4 <= Re < 4e5:
        C = 0.027
        m = 0.805
    Nu = C*Re**m*Pr**(1/3)
    return Nu
    
def Zukauskas():
# Zukauskas method
    if 1 <= Re < 40:
        C = 0.75
        m = 0.4
    elif 40 <= Re < 1e3:
        C = 0.51
        m = 0.5
    elif 1e3 <= Re < 2e5:
        C = 0.26
        m = 0.6
    elif 2e5 <= Re < 1e6:
        C = 0.076
        m = 0.7
    if Pr <= 10:
        n = 0.37
    elif Pr >= 10:
        n = 0.36
    Nu = C*Re**m*Pr**n*(Pr/Prs)**(1/4)
    return Nu

def Churchill_Burnstein():
# Churchill and Burnstein method
    Nu = 0.3+(0.62*Re**(1/2)*Pr**(1/3))/(1+(0.4/Pr)**(2/3))**(1/4)*(1+(Re/282000)**(5/8))**(4/5)
    return Nu

# Free convection correlations
def Morgan():
# Morgan method
    Ra = 9.81*(1/T_f)*(T_s-T_inf)*D**3/(visc**2)*Prf
    if 1e-10 <= Ra < 1e-2:
        C = 0.675
        n = 0.058
    elif 1e-2 <= Ra < 1e2:
        C = 1.02
        n = 0.148
    elif 1e2 <= Ra < 1e4:
        C = 0.850
        n = 0.188
    elif 1e4 <= Ra < 1e7:
        C = 0.480
        n = 0.250
    elif 1e7 <= Ra < 1e12:
        C = 0.125
        n = 0.333
    Nu = C*Ra**n
    return Nu

def Churchill_Chu():
# Churchill and Chu method
    Ra = 9.81*(1/T_f)*(T_s-T_inf)*D**3/(visc**2)*Prf
    Nu = (0.60 + 0.387*Ra**(1/6)/(1+(0.559/Prf)**(9/16))**(8/27))**2
    return Nu

# Calculate Nu and h_bar from all applicable methods 
def print_Nu():
# return Nu for all methods applicable to conditions
    if condition == 'forced':  
        if Re < 0.4 or 4e5 < Re <= 1e6:
            print "Nu from Zukauskas: %0.3f" %Zukauskas()
            print "Nu from Churchill and B: %0.3f" %Churchill_Burnstein()
            return
        elif Re < 1 or Re > 1e6:
            print "Nu from Churchhill and Burnstein: %0.3f" %Churchill_Burnstein()
            return
        elif Re >= 0.4:
            print "Nu from Zukauskas: %0.3f" %Zukauskas()
            print "Nu from Churchill and Burnstein: %0.3f" %Churchill_Burnstein()
            print "Nu from Hilpert: %0.3f" %Hilpert()
            return
    elif condition == 'free':
        print "Nu from Morgan: %0.3f" %Morgan()
        print "Nu from Churchill and Chu: %0.3f" %Churchill_Chu()
        return
       
def print_h_bar():
# return h_bar from all methods applicable 
    if condition == 'forced':    
        if Re < 0.4 or 4e5 < Re <= 1e6:
            hZ = Zukauskas()*(k/D)
            hC = Churchill_Burnstein()*(k/D)
            print "h_bar from Zukauskas: %0.3f" %hZ
            print "h_bar from Churchill and B: %0.3f" %hC
            return
        elif Re < 1 or Re > 1e6:
            hC = Churchill_Burnstein()*(k/D)
            print "h_bar from Zukauskas: %0.3f" %hC
            return
        elif Re >= 0.4:
            hZ = Zukauskas()*(k/D)
            hC = Churchill_Burnstein()*(k/D)
            hH = Hilpert()*(k/D)
            print "h_bar from Zukauskas: %0.3f" %hZ
            print "h_bar from Churchill and Burnstein: %0.3f" %hC
            print "h_bar from Hilpert: %0.3f" %hH
            return
    elif condition == 'free':
        hM = Morgan()*(k/D)
        hC = Churchill_Chu()*(k/D)
        print "h_bar from Morgan: %0.3f" %hM
        print "h_bar from Churchill and Chu: %0.3f" %hM
        return
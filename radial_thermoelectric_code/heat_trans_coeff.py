# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 15:31:24 2016
Calculating heat transfer coefficient for radial thermoelectric project.
author: Michael Orrill
"""
#%% import libraries
from __future__ import division
import numpy as np
import os.path
import itertools
#%% Read in data for air and water
def read_data(file):
    if os.path.isfile(file):
        infile = open(file)
        volt = np.zeros((101,1))
        current = np.zeros((101,1))
        f = csv.reader(infile)
        count = 0
        for row in itertools.islice(f,7,108):
            volt[count] = row[1]
            current[count] = row[13]
            count = count + 1
        infile.close()
        # caclulate slope of data set
        dI = np.gradient(current[:,0])
        slopeV = np.gradient(volt[:,0],dI)
        Res = np.mean(slopeV)
        return volt, current, slopeV, Res
    else:
        #print "File: "+ file[59:] +" not found in directory"
        volt = []
        current = []
        slopeV = []
        Res = []
        return volt, current, slopeV, Res
#%% Input flow conditions and fluid properties

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


#%% Create functions for each correlation
def Hilpert(Re,Pr):
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
    
def Zukauskas(Re,Pr,Prs):
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

def Churchill_Burnstein(Re,Prf):
    Nu = 0.3+(0.62*Re**(1/2)*Pr**(1/3))/(1+(0.4/Pr)**(2/3))**(1/4)*(1+(Re/282000)**(5/8))**(4/5)
    return Nu
# return Nu for all methods applicable to conditions   
def Nu():
    if Re < 0.4:
        print "Nu from Zukauskas: %0.3f" %Zukauskas(Re,Pr,Prs)
        print "Nu from Churchill and B: %0.3f" %Churchill_Burnstein(Re,Prf)
    elif Re < 1:
        print "Nu from Churchhill and Burnstein: %0.3f" %Churchill_Burnstein(Re,Prf)
    elif Re >= 0.4:
        print "Nu from Zukauskas: %0.3f" %Zukauskas(Re,Pr,Prs)
        print "Nu from Churchill and Burnstein: %0.3f" %Churchill_Burnstein(Re,Prf)
        print "Nu from Hilpert: %0.3f" %Hilpert(Re,Pr)
# return h_bar from all methods applicable        
def get_h_bar():
    if Re < 0.4:
        hZ = Zukauskas(Re,Pr,Prs)*(k/D)
        hC = Churchill_Burnstein(Re,Prf)*(k/D)
        print "h_bar from Zukauskas: %0.3f" %hZ
        print "h_bar from Churchill and B: %0.3f" %hC
    elif Re < 1:
        hC = Churchill_Burnstein(Re,Prf)*(k/D)
        print "h_bar from Zukauskas: %0.3f" %Churchill_Burnstein(Re,Pr)*(k/D)
    elif Re >= 0.4:
        hZ = Zukauskas(Re,Pr,Prs)*(k/D)
        hC = Churchill_Burnstein(Re,Prf)*(k/D)
        hH = Hilpert(Re,Pr)*(k/D)
        print "h_bar from Zukauskas: %0.3f" %hZ
        print "h_bar from Churchill and Burnstein: %0.3f" %hC
        print "h_bar from Hilpert: %0.3f" %hH
#%% 



















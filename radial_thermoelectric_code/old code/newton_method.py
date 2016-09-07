# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 12:16:54 2016

@author: Michael Orrill
"""
#%% solver libraries
from scipy.optimize import fsolve
import math as m
import numpy as np
from __future__ import division
#%%
x1 = 1.0
x2 = 1.0
x3 = 1.0
while True:
    #functions
    f1 = x1**3.0 - 2*x2 -2
    f2 = x1**3.0 - 5*x3**2.0 + 7
    f3 = x2*x3**2.0 - 1
    f = np.array(([f1],[f2],[f3]))      
    # Jacobian
    J = np.array(([3*x1**2.0, -2.0, 0], [3*x1**2.0, 0, -10*x3], [0, x3**2.0, 2*x2*x3]))
    
    [x1_next,x2_next,x3_next] = np.array(([x1],[x2],[x3])) - np.matmul(np.linalg.inv(J),f)
    t = abs(x1_next[0]-x1)+abs(x2_next[0]-x2)+abs(x3_next[0]-x3)
    print t
    if t < 0.0001:
        print "Using Newtons method:\nx1 = %0.3f\nx2 = %0.3f\nx3 = %0.3f" %(x1,x2,x3)
        break
    x1,x2,x3 = x1_next[0],x2_next[0],x3_next[0]
def F(p):
    x,y,z = p
    return ((x**3.0 - 2*y -2, x**3.0 - 5*z**2.0 + 7, y*z**2.0 - 1))

[x,y,z] = fsolve(F,(1,1,1))
print "Using fsolve:\nx1 = %0.3f\nx2 = %0.3f\nx3 = %0.3f" %(x,y,z)
#%% trying newtons method with temperature dist.
# geometry
rwi = 0.05          # inner pipe radius [m]
rtei = 0.055         # outer pipe/inner TE radius [m]
rteo = 0.1          # outer TE radius [m]
tn = 0.4            # n-type thickness [m]
tp = 0.4            # p-type thickness [m]
ts = abs(1.0-tn-tp) # spacer thickness [m]
tw = tn+ts+tp       # pipe thickness [m]
# temperature and flow conditions
Th = 300.0          # hot stream tempertaure [K]
Tc = 100.0          # cold stream temperature [K]
hh = 53.216            # hot side convection coefficient [W/m^2.K]
hc = 67.83         # cold side convection coefficien [W/m^2.K]
# material properties
kw = 1.0            # pipe thermal conductivity [W/m.K]
kn = 1.0            # n-type material thermal conductivity [W/m.K]
kp = 1.0            # p-type material thermal conductivity [W/m.K]
ks = 1.0            # spacer/insulator thermal conductivity [W/m.K]
Sn = -1.0            # n-type Seebeck coefficient [µV/K] 
Sp = 1.0            # p-type Seebeck coefficient [µV/K] 
sigma_n = 1.0       # n-type electrical conductivity []
sigma_p = 1.0       # p-type electrical conductivity []

# thermal resistances:
Rte = (m.log(rteo/rtei)/(2.0*m.pi))*(1.0/(kn*tn)+1.0/(ks*ts)+1.0/(kp*tp)) # TE/spacer total R
Rw = m.log(rteo/rtei)/(2.0*m.pi*tw*kw)
# electrical resistance
Re = (m.log(rteo/rtei)/(2.0*m.pi))*(1.0/(sigma_n*tn)+1.0/(sigma_p*tp)) # n- and p-type total Re

Tw = 140
T1 = 135
T2 = 130
while True:
    #functions
    f1 = hh*2.0*m.pi*rwi*tw*(Th-Tw) - 2.0*m.pi*kw*tw/m.log(rtei/rwi)*(Tw-T1)
    f2 = 2.0*m.pi*kw*tw/m.log(rtei/rwi)*(Tw-T1) - (T1-T2)/Rte + (Sp-Sn)**2.0*(T1-T2)**2.0/(4.0*Re) - (Sp-Sn)**2.0*(T1-T2)*T1/(2.0*Re)
    f3 = (T1-T2)/Rte + (Sp-Sn)**2.0*(T1-T2)**2.0/(4.0*Re) + (Sp-Sn)**2.0*(T1-T2)*T2/(2.0*Re) - hc*2.0*m.pi*rteo*tw*(T2-Tc)
    f = np.array(([f1],[f2],[f3]))        
    # Jacobian
    J = np.array(([-hh*2.0*m.pi*rwi*tw-2.0*m.pi*kw*tw/m.log(rtei/rwi), 2.0*m.pi*kw*tw/m.log(rtei/rwi), 0], [2.0*m.pi*kw*tw/m.log(rtei/rwi), -2.0*m.pi*kw*tw/m.log(rtei/rwi)-1/Rte+(Sp-Sn)**2*(2.0*T1-T2)/(4.0*Re)-(Sp-Sn)**2*(2.0*T1-T2)/(2.0*Re), 1/Rte+(Sp-Sn)**2*(2.0*T2-T1)/(4.0*Re)-(Sp-Sn)**2*(-T1)/(2.0*Re)], [0, 1/Rte+(Sp-Sn)**2*(2.0*T1-T2)/(4.0*Re)+(Sp-Sn)**2*(T2)/(2.0*Re), -1/Rte+(Sp-Sn)**2*(2.0*T2-T1)/(4.0*Re)+(Sp-Sn)**2*(T1-2.0*T2)/(2.0*Re)-hc*2.0*m.pi*rteo*tw]))
    [Tw_next,T1_next,T2_next] = np.array(([Tw],[T1],[T2])) - np.matmul(np.linalg.inv(J),f)
    t = abs(Tw_next[0]-Tw)+abs(T1_next[0]-T1)+abs(T2_next[0]-T2) 
    print t
    print (Tw_next[0],T1_next[0],T2_next[0])
    if t < 1e-5:
        print "Tw = %0.3f\nT1 = %0.3f\nT2 = %0.3f" %(Tw,T1,T2)
        break
    Tw,T1,T2 = Tw_next[0],T1_next[0],T2_next[0]

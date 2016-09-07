# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 12:44:41 2016

@author: Michael Orrill
"""
#%% solver libraries
from scipy.optimize import fsolve
import math as m
import numpy as np
from __future__ import division
#%% Variebles/parameters
# geometry
rwi = 0.05          # inner pipe radius [m]
rtei = 0.055         # outer pipe/inner TE radius [m]
rteo = 0.1          # outer TE radius [m]
tn = 0.4            # n-type thickness [m]
tp = 0.4            # p-type thickness [m]
ts = abs(1.0-tn-tp) # spacer thickness [m]
tw = tn+ts+tp       # pipe thickness [m]
# temperature and flow conditions
Th = 475.0          # hot stream tempertaure [K]
Tc = 275.0          # cold stream temperature [K]
hh = 53.216            # hot side convection coefficient [W/m^2.K]
hc = 67.83            # cold side convection coefficien [W/m^2.K]
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
ans()

#%% solver 
def ans():
    Tw,T1,T2 = fsolve(f,(100,100,100))
    print "Method 1:\nTw = %0.3f\nT1 = %0.3f\nT2 = %0.3f" %(Tw,T1,T2)
    Tw,T1,T2 = fsolve(f2,(100,100,100))
    print "Method 2:\nTw = %0.3f\nT1 = %0.3f\nT2 = %0.3f" %(Tw,T1,T2)

#%% funtion for equations
def f(p):
    Tw,T1,T2 = p
    return ((hh*2.0*m.pi*rwi*tw)*(Th-Tw) - 2.0*m.pi*kw*tw/m.log(rtei/rwi)*(Tw-T1), 2.0*m.pi*kw*tw/m.log(rtei/rwi)*(Tw-T1) - (T1-T2)/Rte + (Sp-Sn)**2.0*(T1-T2)**2.0/(4.0*Re) - (Sp-Sn)**2.0*(T1-T2)*T1/(2.0*Re), (T1-T2)/Rte + (Sp-Sn)**2.0*(T1-T2)**2.0/(4.0*Re) + (Sp-Sn)**2.0*(T1-T2)*T2/(2.0*Re) - hc*2.0*m.pi*rteo*tw*(T2-Tc))

Tw,T1,T2 = fsolve(f,(1,1,1))
print "Tw = %0.3f\nT1 = %0.3f\nT2 = %0.3f" %(Tw,T1,T2)

# does the same as above a little more clearly
def f2(p):
    Tw = p[0]
    T1 = p[1]
    T2 = p[2]
    F = np.zeros((3))
    F[0] = hh*2.0*m.pi*rwi*tw*(Th-Tw) - 2.0*m.pi*kw*tw/m.log(rtei/rwi)*(Tw-T1)
    F[1] = 2.0*m.pi*kw*tw/m.log(rtei/rwi)*(Tw-T1) - (T1-T2)/Rte + (Sp-Sn)**2.0*(T1-T2)**2.0/(4.0*Re) - (Sp-Sn)**2.0*(T1-T2)*T1/(2.0*Re)
    F[2] = (T1-T2)/Rte + (Sp-Sn)**2.0*(T1-T2)**2.0/(4.0*Re) + (Sp-Sn)**2.0*(T1-T2)*T2/(2.0*Re) - hc*2.0*m.pi*rteo*tw*(T2-Tc)
    return F
Tw,T1,T2 = fsolve(f2,(1,1,1))
print (Tw,T1,T2)

# manual Newton's method:
def newton(Tw0,T10,T20):
    Tw = Tw0
    T1 = T10
    T2 = T20
    while True:
        #functions
        f1 = hh*2.0*m.pi*rwi*tw*(Th-Tw) - 2.0*m.pi*kw*tw/m.log(rtei/rwi)*(Tw-T1)
        f2 = 2.0*m.pi*kw*tw/m.log(rtei/rwi)*(Tw-T1) - (T1-T2)/Rte + (Sp-Sn)**2.0*(T1-T2)**2.0/(4.0*Re) - (Sp-Sn)**2.0*(T1-T2)*T1/(2.0*Re)
        f3 = (T1-T2)/Rte + (Sp-Sn)**2.0*(T1-T2)**2.0/(4.0*Re) + (Sp-Sn)**2.0*(T1-T2)*T2/(2.0*Re) - hc*2.0*m.pi*rteo*tw*(T2-Tc)
        f = np.array(([f1],[f2],[f3]))        
        # Jacobian
        J = np.array(([-hh*2.0*m.pi*rwi*tw-2.0*m.pi*kw*tw/m.log(rtei/rwi), 2.0*m.pi*kw*tw/m.log(rtei/rwi), 0], [2.0*m.pi*kw*tw/m.log(rtei/rwi), -2.0*m.pi*kw*tw/m.log(rtei/rwi)-1/Rte+(Sp-Sn)**2*(2.0*T1-T2)/(4.0*Re)-(Sp-Sn)**2*(2.0*T1-T2)/(2.0*Re), 1/Rte+(Sp-Sn)**2*(2.0*T2-T1)/(4.0*Re)-(Sp-Sn)**2*(-T1)/(2.0*Re)], [0, 1/Rte+(Sp-Sn)**2*(2.0*T1-T2)/(4.0*Re)+(Sp-Sn)**2*(T2)/(2.0*Re), -1/Rte+(Sp-Sn)**2*(2.0*T2-T1)/(4.0*Re)+(Sp-Sn)**2*(T1-2.0*T2)/(2.0*Re)-hc*2.0*m.pi*rteo*tw]))
        
        [Tw_next,T1_next,T2_next] = np.array((Tw,T1,T2)) - np.matmul(np.linalg.inv(J),f)
        t = abs(Tw_next[0]-Tw)+abs(T1_next[0]-T1)+abs(T2_next[0]-T2) 
        print t
        if t < 0.0001:
            print "Tw = %0.3f\nT1 = %0.3f\nT2 = %0.3f" %(Tw,T1,T2)
            break
        Tw,T1,T2 = Tw_next[0],T1_next[0],T2_next[0]
    return Tw,T1,T2
    
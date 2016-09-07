
import numpy as np
import matplotlib.pyplot as plt
from __future__ import division
#%%
start = np.array([-2,-3,-4,-5,-6,-7,-8,-9,-10])
end = -1*start
x1 = np.linspace(start,end,num=100,endpoint=True)
#%%
x1 = np.arange(0,2*np.pi, 0.1)
y1 = 19*np.sin(x1)**3
y2 = 17*np.cos(x1)-8*np.cos(2*x1)-5*np.cos(3*x1)-np.cos(4*x1)
    
plt.figure(1)
plt.plot(y1,y2,'r')
#%%
#x1 = np.linspace(-2,2,num=100,endpoint=True)
u = np.arange(0,19)

for i in range(len(u)):
    x1 = np.arange(0,2*np.pi, 0.1)
    y1 = u[i]*(16/13)*np.sin(x1)**3
    y2 = u[i]*np.cos(x1)-5*np.cos(2*x1)-2*np.cos(3*x1)-np.cos(4*x1)
        
    plt.figure(1,figsize = (10,10))
    plt.plot(y1,y2,'r')
    plt.grid()
plt.show()

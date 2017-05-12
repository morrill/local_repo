# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 11:11:16 2017

@author: Michael Orrill
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import ScalarFormatter
rcParams['font.size'] = 25
#%%
c = np.loadtxt('C4.1_no_sonicate.txt')
fig = plt.figure(figsize=(20,10))
ax = fig.add_subplot(111)
ax.semilogx(c[:,0],c[:,1],color='r',linewidth=3)
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.set_title('Particle size distribution of HCNS in ethylene glycol')
ax.set_ylabel('Intensity [%]')
ax.set_xlabel('Diameter [nm]')

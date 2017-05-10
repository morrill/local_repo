# -*- coding: utf-8 -*-
"""
Created on Wed May 10 11:23:42 2017

@author: Michael Orrill
"""

# this if for plotting the values of Z = 1/Oh in a We/Re space
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib import rcParams
rcParams['font.size'] = 14
#%% Create arrays
re = np.arange(0.1,1000)
we = np.arange(1,1001)
lowerbound = 1  # set as the lower Z boundary line to be shown in plot
upperbound = 10 # set as the upper Z boundary line to be shown in plot

# actual functions to be plotted
low_line = re**2/lowerbound**2
high_line = re**2/upperbound**2

#%% Create plot
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
ax.loglog(re, low_line, color='r', linestyle='--', linewidth=2)
ax.loglog(re, high_line, color='r', linestyle='--', linewidth=2)
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.yaxis.set_major_formatter(ScalarFormatter())
plt.axhline(4, color='k', linestyle='--')
plt.xlim(0.1,1000)
plt.ylim(1,1500)
plt.xlabel('Reynolds number', fontweight='bold')
plt.ylabel('Weber number', fontweight='bold')

plt.fill_between(low_line,[0,1000], color='b')

#ax.xaxis.set_major_formatter(ScalarFormatter())
#ax.yaxis.set_major_formatter(ScalarFormatter())

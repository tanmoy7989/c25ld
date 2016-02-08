#!/usr/bin/env python
import numpy as np
from scipy import integrate
import os, sys
import matplotlib.pyplot as plt

# user -input
fftype = sys.argv[1]
data_dir = sys.argv[2]
datalog_prefix = sys.argv[3]

# GLOBALS
#LDCuts = np.arange(4.5, 10.0 , 0.1)

def getAnalyticalSrel(ax):
    deriv_data = np.loadtxt(os.path.join(data_dir, '%s_dSreldrc.dat' % datalog_prefix))
    Nvals = len(deriv_data)
    srel = np.zeros([Nvals], np.float64)
    
    #srel[1:] = integrate.cumtrapz(deriv_data[:,1], deriv_data[:,0])
    
    current_srel = 0.0
    ind = 0
    for this_srel in deriv_data[:-1,1]:
        current_srel += (this_srel * 0.1)
        srel[ind+1] = current_srel
        ind += 1
   
    ax.plot(deriv_data[:,0], srel[0:], linestyle = 'solid', color = 'blue', linewidth = 3,
            marker = 'o', markerfacecolor = 'blue', markeredgecolor = 'black',
            label = 'analytical')
        
    return srel
    

def getBennettSrel(ax_b, ax_fep):
    data_bennett = np.loadtxt(os.path.join(data_dir, '%s_Srel_Bennett.dat' % datalog_prefix))
    srel_bennett = np.zeros(len(data_bennett), np.float64)
    
    data_fep = np.loadtxt(os.path.join(data_dir, '%s_Srel_Fep.dat' % datalog_prefix))
    srel_fep = np.zeros(len(data_fep), np.float64)
    
    current_srel = 0 
    for i, this_srel in enumerate(data_bennett[:,1]):
        current_srel += this_srel
        srel_bennett[i] = current_srel
        
    current_srel = 0 
    for i, this_srel in enumerate(data_fep[:,1]):
        current_srel += this_srel
        srel_fep[i] = current_srel
    
    
    ax_b.plot(data_bennett[:,0], srel_bennett[0:], linestyle = 'solid', color = 'green', linewidth = 3,
             marker = 'o', markerfacecolor = 'green', markeredgecolor = 'black',
             label = 'Bennett') 

    ax_fep.plot(data_fep[:,0], srel_fep[0:], linestyle = 'solid', color = 'green', linewidth = 3,
                marker = 'o', markerfacecolor = 'green', markeredgecolor = 'black',
                label = 'FreeEnePerturb')  
            
    #ind = np.argmin(srel_bennett)
    #print ind, LDCuts[ind]
    return srel_bennett, srel_fep
    

fig = plt.figure(figsize = (10,4), facecolor = 'w', edgecolor = 'w')
ax1 = fig.add_subplot(1,3,1)
ax2 = fig.add_subplot(1,3,2)
ax3 = fig.add_subplot(1,3,3)

AnalyticalSrel = getAnalyticalSrel(ax1)
BennettSrel, FepSrel = getBennettSrel(ax2, ax3); 
                   
## Design               
for ax in [ax1, ax2, ax3]:
    ax.legend(loc = 'best', prop = {'size': 'small'})
    ax.set_xlabel(r'$r_c (\AA)$', fontsize = 'large')
    ax.set_ylabel(r'$S_{rel}$' , fontsize = 'large')

plt.subplots_adjust(left = 0.1, bottom = 0.2, wspace = 0.5)
plt.show()

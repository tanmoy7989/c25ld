#!/usr/bin/env python

import os, sys, pickle
import numpy as np
import matplotlib.pyplot as plt

pickleName = sys.argv[1]
data = pickle.load(open(pickleName, 'r'))

ravg = data[0]
NCuts = len(ravg)
sigmasq = np.zeros(NCuts, np.float64)
sigmasq_bulk = np.zeros(NCuts, np.float64)

kB = 1.3806e-23
Temp = 298.0
kappa = 4.5227e-10
rho = (1500/(35.50635e-10)**3.)
prefac = (kB * Temp) * kappa * rho**2

ravg = np.array(ravg)
ravg *= 0.9
LD = data[1]

for j in range(NCuts):
	this_LD = LD[:,j]
	sigmasq[j] = np.var(this_LD)
	sigmasq_bulk[j]  = prefac * (4*np.pi/3.) * ravg[j]**3.
	
print sigmasq_bulk
plt.plot(ravg, sigmasq, linestyle = 'solid', linewidth = 2, color = 'blue',
		 marker = 'o', markersize = 8, markeredgecolor = 'blue', label = 'simulation')
		 
#plt.plot(ravg, sigmasq_bulk*1e-32, linestyle = 'solid', linewidth = 2, color = 'red',
#		 marker = 'o', markersize = 8, markeredgecolor = 'red', label = 'theory')


plt.xlabel(r'$ \mathrm{Cutoff } \langle r_c \rangle (A^o)$', fontsize = 15)
plt.ylabel(r'$\sigma^2(\rho_{LD})$', fontsize = 15)
plt.title('Local Density fluctuations in SPC/E water')

plt.show()

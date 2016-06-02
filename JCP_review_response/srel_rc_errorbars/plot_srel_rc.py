#!/usr/bin/env python

import sys, copy
import numpy as np
import matplotlib.pyplot as plt

system = sys.argv[1]

LDCuts = np.arange(4.5, 10.0, 0.1)
data = np.loadtxt('%s_wca_Srel_Bennett_N4.dat' % system)
NBlocks = data.shape[1]

# mean value to be plotted
srel = np.zeros([len(LDCuts), NBlocks], np.float64)
curr_srel = np.zeros(NBlocks)
for i, x in enumerate(LDCuts):
    curr_srel += data[i,:]
    srel[i,:] = curr_srel

srel_mean = np.mean(srel, axis = 1)
err = np.std(data, axis = 1, ddof = 1)
err_final = np.zeros(len(LDCuts), np.float64)
for i, x in enumerate(LDCuts[1:-1]): err_final[i] = 0.5 * (err[i-1] + err[i+1])

plt.figure()
plt.errorbar(LDCuts, srel_mean, yerr = err_final)
plt.xlabel(r'$r_c$', fontsize = 14)
plt.ylabel(r'$S_{rel}$', fontsize = 14)
plt.title('Block errors')

plt.show()

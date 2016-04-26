#!/usr/bin/env python

import matplotlib.pyplot as plt
import os, sys, pickle
import numpy as np
from matplotlib.ticker import MaxNLocator

target_dir = sys.argv[1]


# Look in write_LD_sensitivity_jobs.py
lbls = ['wca']
#lbls = ['spce']

nrows = 3
ncols = 5
NCuts = 14

r = np.zeros([NCuts,2], np.float64)
r_np = np.zeros([NCuts,2], np.float64)

for i, lbl in enumerate(lbls):
	fig = plt.figure(figsize = (15,3), facecolor = 'white', edgecolor = 'white')
	pickleName = os.path.abspath(os.path.join(target_dir, 'firstshell_%s.pickle' % lbl))
	data = pickle.load(open(pickleName, 'r'))
	LDCuts = data[0]
	LD = data[1]
	NW = data[2]
	splot_count = 1
	for j in range(NCuts):
		this_LD = LD[:,j]
		this_NW = NW
		#ind = np.argsort(this_LD);
		#this_LD = np.sort(this_LD)
		#this_NW = NW[ind]
		
		# calculating corrrelation statistics
		Nobs = len(this_NW)
		mu_NW = np.mean(this_NW); mu_LD = np.mean(this_LD)
		sigma_NW = np.std(this_NW); sigma_LD  = np.std(this_LD)
		for k in range(Nobs):
			r[j,i] += this_LD[k] * this_NW[k]
		r[j,i] -= Nobs * mu_NW * mu_LD
		r[j,i] /= (Nobs - 1) * sigma_NW * sigma_LD	
		
		#if (j==0 or j==2 or j==4 or j==8):		
		if (1):
			ax = plt.subplot(nrows, ncols, splot_count)		
			splot_count +=1
			ax.plot(this_LD,this_NW, linestyle = 'None', linewidth = 2, 
				 	marker = 'x', markersize = 5, markeredgecolor = 'red',
				 	color = 'red', label = 'LD Cut = %g' % LDCuts[j])
		

			ax.xaxis.set_major_locator(MaxNLocator(5, prune = 'both'))
			ax.legend(loc = 'best', prop = {'size': 15}, fancybox = True)
			if j == 0:
				ax.set_ylabel('first shell waters', fontsize = 15)
				ax.set_title('forcefield : %s' % lbl, fontsize = 10)
			else:
				ax.set_ylabel('')
				ax.set_yticklabels('')
	plt.subplots_adjust(wspace = 0.001)		

# save data to pickle
pickleName = os.path.abspath(os.path.join(target_dir, 'corrcoeff.pickle'))
pickle.dump((LDCuts, r), open(pickleName, 'w'))

# plot r^2
plt.figure(figsize = (4,4))
clrs = ['red', 'blue']
maxcorr_cut = 6.5
for i,lbl in enumerate(lbls):
	plt.plot(LDCuts, r[:,i]**2, linestyle = 'solid', marker = 'o',
			color = clrs[i], markerfacecolor = clrs[i], 
			linewidth = 3, markersize = 8, label = lbl)
	plt.plot(maxcorr_cut*np.ones(10), np.linspace(0, 0.98 * (r**2.).max(),10), linestyle = 'dashed', color = 'black', linewidth = 2) 
	
	print LDCuts, r[:,i]**2
	R_max = np.max(r[:,i]**2)
	maxind = np.where(r[:,i]**2 == R_max)[0]
	cut_max = LDCuts[maxind]
	print lbl, "   R_max, Cut_max",  R_max, cut_max
	plt.xlabel('local density cut' + r'$<r_c> A^o$', fontsize = 15)
	plt.ylabel(r'$R^2$', fontsize = 15)
	

plt.title('Squared correlation between first shell waters (fsw) and local density (ld)')
plt.legend()
																									
plt.show()
		
		

#/usr/bin/env python

import numpy as np
import os, sys, pickle
import matplotlib.pyplot as plt

prefix = 'c25'
lbls = ['lj', 'wca']
rdfs = ['MM', 'WW', 'MW']
ind = '22'

target_dir = sys.argv[1]

prefix  = 'c25_poly'
rdfs = ['MM']; ind = '12'

def getCuts(r, rdf):
	Npoints = len(rdf)
	cuts = np.array([])
	for i in range(1,Npoints-1):
		if rdf[i] < rdf[i-1] and rdf[i] < rdf[i+1]:
			cuts = np.append(cuts, r[i])
	
	## Null-check 04/30/15
	if not cuts:
		cuts = np.append(cuts, r[0])
	return cuts


for i, lbl in enumerate(lbls):
	pickleName = os.path.abspath(os.path.join(target_dir, '%s_%s_rdf.pickle' % (prefix, lbl)))
	data = pickle.load(open(pickleName, 'r'))
	x = data['bin_centers']

	plt.figure(figsize = (4,4), facecolor = 'white', edgecolor = 'white')
	
	for j, rdf_type in enumerate(rdfs):	
		key = rdf_type
		y = data[key]['g']
		fscut = getCuts(x,y)[0]
		print rdf_type , ':', fscut
		ax = plt.subplot(ind + str(j+1))
		ax.plot(x,y,linestyle = 'solid', linewidth = 3, color = 'red', 
				marker = 'o', markerfacecolor = 'black', markeredgecolor = 'black', 
				label = r'$g(r) \mathrm{ %s, %s}$' % (rdf_type, lbl))
	
		ax.plot(fscut*np.ones(10), np.linspace(0, 0.98 * y.max(),10), linestyle = 'dashed', linewidth = 2, color = 'black')
		
		ax.set_xlabel(r'$r$', fontsize = 15)
		ax.set_ylabel(r'$g(r)$', fontsize = 15)
		ax.legend(loc = 'best')
		plt.subplots_adjust(bottom = 0.2, left = 0.2)
	plt.title(lbl)
plt.show()

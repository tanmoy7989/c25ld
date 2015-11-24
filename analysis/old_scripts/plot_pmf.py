#!/usr/bin/env python

import os, sys, pickle
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

target_dir = os.path.abspath(sys.argv[1])
fftype = sys.argv[2]

#MATPLOTLIB rc params
matplotlib.rc('figure', facecolor = 'w', edgecolor = 'w', figsize = (4,4))
matplotlib.rc('lines', linewidth = 3, marker = 'None', markersize = 8)

methods = ['pymbar_wham', 'grossfield_wham']
clrs = ['red', 'blue', 'green', 'black']
showErr = True

# Unit conversions
kB = 1.38065e-23 #J/Kelvin
N_A = 6.023e23 #Avogadro's number
Temp = 298.0 #Kelvin
conv_factor = (kB *N_A *1e-3) * Temp
ref_global = 0.

for method in methods:
	fig = plt.figure()
	ax = plt.subplot('111')
	prefix_list = ['AA_'+fftype, 
				   'CG_'+fftype+'_SP', 
			  	   'CG_'+fftype+'_SPLD', 
			  	   'CG_'+fftype+'_LD']
	for i, prefix in enumerate(prefix_list):
		pickleName = prefix + '_' + method + '_pmf.pickle'
		pickleName = os.path.join(target_dir, pickleName)
		try:
			data = pickle.load(open(pickleName, 'r'))
		except IOError:
			print 'WHAM pickle for %s not found' % prefix
			continue
		x = data['bin_centers'] - data['bin_centers'].min(); x/= 10
		y = data['pmf'] * conv_factor # convert pmf to kJ/mol-K
		err = data['std'] * conv_factor
		
		## fix a reference, 05/26/2015, TS
		if prefix.startswith('AA'):
			displacement = 0.
			ref_global = y.min()
		elif prefix.startswith('CG'):
			displacement = abs(y.min() - ref_global)
		
		y -= displacement
		if showErr:
			ax.errorbar(x,y,yerr = err, color = clrs[i], label = prefix + '_'+ method.split('_')[0])
		else:
			ax.plot(x,y,color = clrs[i], label = prefix + '_'+ method.split('_')[0])
					
	ax.set_xlabel(r'$Rg - Rg^{min} \mathrm{ (nm)}$', fontsize = 15)
	ax.set_ylabel(r'$W(Rg) (kJ/mol)$', fontsize = 15)
	leg = ax.legend(loc = 'best', prop = {'size': 15})
	leg.get_frame().set_alpha(0.2)

	ax.set_xlim([0,1])	
	# Plot-design
	plt.subplots_adjust(left = 0.2, bottom = 0.2)
	if showErr:
		MaxErr = np.max(err/y) 
		print 'Relative Errors for ', method
		print (abs(err)/abs(y)) 
	
	# Saving to disk
	figname = os.path.join(target_dir, 'pmf_%s_%s.svg' % (fftype, method))
	plt.savefig(figname, bbox_inches = 'tight', dpi = 300)

plt.show()
	
	

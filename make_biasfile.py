#!/usr/bin/env python

import os
import numpy as np

def make_Biasfile(targetdir, filename):
	
	# Uniform Biasing
	
	def uniform_bias(Ncenters, Rg_min, Rg_max, Fconst, fname):
		Rg = np.linspace(Rg_min, Rg_max, Ncenters)
		fconst = Fconst*np.ones([Ncenters])
		f = open(fname, 'w')
		for i in range(Ncenters):
			f.write('%f\t%f\n' % (Rg[i], fconst[i]))
		f.close()
	
	Rg_min = 2.0
	Rg_max = 18.0
	Ncenters = 25
	
	# Large Bias chosen
	Fconst = 32.00 #kcal/mol/A^2
	# note that this is the 'k' for u_bias = (k/2)*(x-x0)^2
	
	fname = os.path.join(targetdir, filename)
	uniform_bias(Ncenters, Rg_min, Rg_max, Fconst, fname)
	return fname


import sys
targetdir = sys.argv[1]
filename = sys.argv[2]
make_Biasfile(targetdir, filename)	

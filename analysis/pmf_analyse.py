#usr/bin/env python
import numpy as np
import os, pickle, sys

sys.path.append(os.path.expanduser('~/c25ld/analysis'))
from utils import Parsers
from mywhamlib import Whamlib

import matplotlib.pyplot as plt

# User-inputs
trajdir = os.path.abspath(sys.argv[1])
paramfile = os.path.abspath(sys.argv[2])
trjsrc = sys.argv[3]
trjprefix = sys.argv[4]
prefix = sys.argv[5]
target_dir = sys.argv[6]

methods = ['pymbar_wham', 'grossfield_wham']
nbins = 50
block_dict = {'lmp': 1, 'sim': 4}
nblocks = block_dict[trjsrc]
doSubSample = True

for method in methods:
	PMFpickle = prefix + '_' + method + '_pmf.pickle'
	PMFpickle = os.path.join(target_dir, PMFpickle)
	
	
	# check if data already present
	if os.path.isfile(PMFpickle):
			print 'PMF pickle already present. Skipping computation...'
			continue

	print 'Computing umbrella sampled pmf for %s using %s'% (prefix, method)
	# Set parameters for parsers
	nstages = len(np.loadtxt(paramfile)[:,0])
	Parsers.nstages = nstages
	Parsers.Trajdir = trajdir
	Parsers.Prefix = trjprefix
	Parsers.trjsrc = trjsrc
	fp = Parsers.parse_PElog
	fu = Parsers.parse_Ulog

	# Bypass whamlib's internal parsers and parse out all data at once
	## 04/16/15: Global pickle loading enabled
	
	global_whamdatapickle = prefix + '_global_whamdata.pickle'
	if os.path.isfile(global_whamdatapickle):
		print 'Loading from global data pickle'
		(E_pe, E_u, Rg, nframes_global) = pickle.load(open(global_whamdatapickle, 'r'))
	else:
		E_pe = fp()
		(E_u, Rg) = fu()	
		nframes_global = E_pe.shape[1]
		pickle.dump((E_pe,E_u,Rg,nframes_global), open(global_whamdatapickle, 'w'))
		
	# Calling WHAM routines
	nframes_per_block = int(nframes_global/nblocks)
	pmf_block = np.zeros([nbins,nblocks], np.float64)
	bin_centers = np.zeros(nbins, np.float64) 
	count = 0
	
	# block iteration
	for block in range(nblocks):
		if nblocks == 1:	
			whamPrefix = prefix
		else:
			whamPrefix = prefix + '_block' + str(block)
		print "Running WHAM on ", whamPrefix
		
		Step = {'nstages': nstages, 'nframes': nframes_per_block}
		wham = Whamlib.WHAM(Prefix = whamPrefix, TrajPrefix = trjprefix, 
							trajdir = trajdir, biasfile = paramfile, 
							Step = Step, nbins = nbins)
		wham.Method = method
	
		# make whamlib read data block by block	
		whamhistpickle = whamPrefix + '_wham_hist.pickle'
		if not os.path.isfile(whamhistpickle):
			wham.E_pe = E_pe[:,count:count+nframes_per_block]
			wham.E_u = E_u[:,count:count+nframes_per_block]
			wham.Rg = Rg[:,count:count+nframes_per_block]
		
		# disable internal parsing of whamlib
		wham.useParser = False
		
		# Subsample data
		## Note that if data is loaded back from a histpickle,
		## frame reduction due to subsampling is already taken into account
		if doSubSample: 
			wham.doSubSample = doSubSample
			wham.freq = 10
		wham.genHist()
		
		#check overlap
		fig = wham.plot_overlap()
		
		# compute pmf
		x = wham.compute_pmf()
		if block == 0:	bin_centers = x[0]
		pmf_block[:,block] = x[1][:]
		count += nframes_per_block
		
	# Error analysis
	pmf = np.sum(pmf_block, axis = 1)/nblocks
	std = np.std(pmf_block,axis=1)
	
	
	# Dumping data to a pickle
	pmf_data = {'bin_centers': bin_centers, 'pmf': pmf, 'std': std,'pmf_block': pmf_block}
	pickle.dump(pmf_data, open(PMFpickle,'w'))

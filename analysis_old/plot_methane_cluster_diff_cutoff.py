#!/usr/bin/env python

import os, sys, pickle
import numpy as np
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt

AAdir = os.path.abspath(sys.argv[1])
CGdir = os.path.abspath(sys.argv[2])
savedir= os.path.abspath(sys.argv[2])

fftypes = ['wca', 'lj']
dirtypes = [os.path.join(CGdir, 'rc'), os.path.join(CGdir, 'rc-1'), os.path.join(CGdir, 'rc+1')]
cgtypes = ['AA', 'CG_SP', 'CG_SPLD']

nrows = 1 ; ncols = len(fftypes)
clrs = {'AA': 'red', 'CG_SP': 'blue', 'CG_SPLD': 'green'}
linest = {0: 'solid', 1: ':', 2:'--'}
dirlabels = {0: 'rc', 1: 'rc-1', 2: 'rc+1'}

axlist = []
fig = plt.figure(figsize = (8,4), facecolor = 'w', edgecolor = 'w')

for i, fftype in enumerate(fftypes):
	ax = fig.add_subplot(nrows, ncols, i+1)
	axlist.append(ax)
	for j, dirtype in enumerate(dirtypes):
		for k, cgtype in enumerate(cgtypes):
			target_dir = ''
			if cgtype == 'AA':
				target_dir = AAdir
				prefix = 'AA_%s_cluster' % fftype
			else:
				target_dir = os.path.join(CGdir, dirtype)
				prefix = 'CG_%s_%s_cluster' % (fftype, cgtype.split('_')[1])
			pickleName = os.path.join(target_dir, 'methane_%s.pickle' % prefix)
			#print pickleName; raw_input()
			try:
				data = pickle.load(open(pickleName, 'r'))
			except IOError:
				print 'Data pickle for %s not present' % prefix
				continue
			x = data[0]
			y = data[1]; y /= np.sum(y)
			
			## Label cgtypes in 1 panel and dirtypes in 2nd panel
			if fftype == 'lj' and cgtype == 'CG_SP':
				lbl = dirlabels[j]
			elif fftype =='wca' and j==0:
				lbl = cgtype
			else:
				lbl = ''
				
			if cgtype == 'AA':
				ax.plot(x,y, linestyle = 'none', color = clrs[cgtype], marker = 'o', markersize = 10, label = lbl)
			else:
				ax.plot(x,y, color = clrs[cgtype], linestyle = linest[j], linewidth = 3, label = lbl)
			ax.set_title(fftype.upper())
		
	
## Design
for i, ax in enumerate(axlist):
	leg = ax.legend(loc = 'best', prop = {'size': 'large'})
	leg.get_frame().set_alpha(0.3)
	ax.set_xlabel('cluster size', fontsize = 'large')
	if i == 0: ax.set_ylabel('distribution', fontsize = 'large')
	ax.yaxis.set_major_locator(MaxNLocator(nbins = 8, prune = 'both'))
	ax.xaxis.set_major_locator(MaxNLocator(nbins = 8, prune = 'both'))
	
plt.subplots_adjust(left = 0.2, bottom = 0.2, hspace = 0.2)
figname = os.path.join(savedir, 'methane_cluster_diff_cutoffs.svg')
fig.savefig(figname, dpi = 300, bbox_inches = 'tight')

plt.show()

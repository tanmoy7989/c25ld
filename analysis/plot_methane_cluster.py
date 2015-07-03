#!/usr/bin/env python

import os, sys, pickle
import numpy as np
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt

target_dir = sys.argv[1]
savedir= sys.argv[2]
fftypes = ['wca']
nrows = len(fftypes); ncols = 1

clrs = ['red', 'blue', 'green', 'cyan']
axlist = []
fig = plt.figure(figsize = (4,4), facecolor = 'w', edgecolor = 'w')

for ind, fftype in enumerate(fftypes):
	prefix_list = ['AA_%s' % fftype, 'CG_%s_SP' % fftype, 'CG_%s_SPLD' % fftype, 'CG_%s_LD' % fftype]
	lbl_list = ['AA', 'CG_SP', 'CG_SPLD', 'CG_LD']
	for i, prefix in enumerate(prefix_list):
		pickleName = os.path.join(target_dir, 'methane_%s_cluster.pickle' % prefix)
		try:
			data = pickle.load(open(pickleName, 'r'))
		except IOError:
			print 'Data pickle for %s not present' % prefix
			continue
			
		x = data[0]
		y = data[1]; y /= np.sum(y)
		ax = fig.add_subplot(nrows,ncols,ind+1)
		if prefix.startswith('AA'):
			ax.plot(x,y,linestyle = 'None', marker = 'o', markersize = 10, color = clrs[i], markeredgecolor = clrs[i], label = lbl_list[i])
		else:
			ax.plot(x,y, linewidth = 3, color = clrs[i], label = lbl_list[i])
		axlist.append(ax)
		
	
## Design
for i, ax in enumerate(axlist):
	if i == 0:
		leg = ax.legend(loc = 'best', prop = {'size': 'x-large'})
		leg.get_frame().set_alpha(0.3)
	
	ax.set_xlabel('cluster size', fontsize = 'xx-large')
	ax.set_ylabel('distribution', fontsize = 'xx-large')
	ax.yaxis.set_major_locator(MaxNLocator(nbins = 8, prune = 'both'))
	
plt.subplots_adjust(left = 0.2, bottom = 0.2, hspace = 0.2)
figname = os.path.join(savedir, 'methane_cluster.svg')
fig.savefig(figname, dpi = 300, bbox_inches = 'tight')

plt.show()

#!/usr/bin/env python

'''
Demonstrate the corrective effect of SPLD over SP
in case of dry and solvated c25 systems. Expect to see that 
said correction is required for the solvated system, but not 
for dry polymer
'''

import os, sys
import numpy as np
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import parse_potential as pp

polydir = os.path.expanduser('~/c25ld/data/cg_ff/polymer')
trajdir = os.path.expanduser('~/c25ld/data/cg_ff/feb15_runs_fsw')

fftype = 'lj'
pAA_label = '%s_AA' % fftype.upper()

polySP_prefix = os.path.join(polydir, 'c25_poly_%s_SP' % fftype)
polySPLD_prefix = os.path.join(polydir, 'c25_poly_%s_SPLD' % fftype)
trajSP_prefix = os.path.join(trajdir, 'c25_%s_SP' % fftype)
trajSPLD_prefix = os.path.join(trajdir, 'c25_%s_SPLD' % fftype)
polyLJ_prefix = os.path.join(polydir, 'c25_poly_%s_LD' % fftype)
trajLJ_prefix = os.path.join(trajdir, 'c25_%s_LD' % fftype)

polySP_logfile, poly_histfile = pp.parseFileNames(polySP_prefix)
polySPLD_logfile, polySPLD_histfile = pp.parseFileNames(polySPLD_prefix)
trajSP_logfile, trajSP_histfile = pp.parseFileNames(trajSP_prefix)
trajSPLD_logfile, trajSPLD_histfile = pp.parseFileNames(trajSPLD_prefix)
polyLJ_logfile, polyLJ_histfile = pp.parseFileNames(polyLJ_prefix)
trajLJ_logfile, trajLJ_histfile = pp.parseFileNames(trajLJ_prefix)

p_polyLJ = pp.parseLog(polyLJ_logfile)
p_trajLJ = pp.parseLog(trajLJ_logfile)
p_polySP = pp.parseLog(polySP_logfile)
p_polySPLD = pp.parseLog(polySPLD_logfile)
p_trajSP = pp.parseLog(trajSP_logfile)
p_trajSPLD = pp.parseLog(trajSPLD_logfile)

fig = plt.figure(figsize = (5,3), facecolor = 'w', edgecolor = 'w')
ax1 = plt.subplot(1,2,1)
ax2 = plt.subplot(1,2,2)

ax1.plot(p_polyLJ[pAA_label][0], p_polyLJ[pAA_label][1], color = 'red', linestyle = 'solid', linewidth = 3, label = 'AA')
ax1.plot(p_polySP['SP'][0], p_polySP['SP'][1], color = 'blue', linestyle = 'solid', linewidth = 3, label = 'dry CG')
ax1.plot(p_trajSP['SP'][0], p_trajSP['SP'][1], color = 'green', linestyle = 'solid', linewidth = 3, label = 'solvated CG')
ax1.set_title('AA vs SP', fontweight = 'bold')

ax2.plot(p_polyLJ[pAA_label][0], p_polyLJ[pAA_label][1], color = 'red', linestyle = 'solid', linewidth = 3, label = 'AA')
ax2.plot(p_polySPLD['SP'][0], p_polySPLD['SP'][1], color = 'blue', linestyle = 'solid', linewidth = 3, label = 'dry CG')
ax2.plot(p_trajSPLD['SP'][0], p_trajSP['SP'][1], color = 'green', linestyle = 'solid', linewidth = 3, label = 'solvated CG')
ax2.set_title('AA vs SPLD', fontweight = 'bold')

for ax in [ax1, ax2]:
	if ax == ax1: 
		leg = ax.legend(loc = 'best', prop = {'size': 12})
		leg.get_frame().set_alpha(0.2)
		ax.set_ylabel(r'$u(r)$' + ' (kcal/mol)', fontsize = 18)
	if ax == ax2:	
		ax.set_yticklabels([])
		
	ax.set_xlabel(r'$r(\AA)$', fontsize = 18)
	ax.xaxis.set_major_locator(MaxNLocator(nbins = 4, prune = 'both'))
	ax.set_ylim([-0.2,0.2])
	
plt.subplots_adjust(left = 0.2, bottom = 0.2, wspace = 0)

plt.show()

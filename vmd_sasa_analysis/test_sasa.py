#!/usr/bin/env python

import numpy as np
import os, sys, pickle
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# Run vmd
def run_vmd(traj, vmd_outprefix, mon_start = 0, mon_end = 25, mon_type = 0):
	if not os.path.isfile(traj):
		print '%s not found\n' % traj
		return

	vmd_params = {'tclscript': 'sasa.tcl', 
	      	      'serial_start': mon_start, 'serial_end': mon_end,
	      	      'start': 0, 'end': -1, 'freq': 100,
	      	      'trajfile': traj,
		      'outprefix': vmd_outprefix,
		      'mon_type' : mon_type}

	vmdstring = 'vmd -dispdev text -e %(tclscript)s -args %(trajfile)s %(outprefix)s %(serial_start)d %(serial_end)d %(start)d %(end)d %(freq)d %(mon_type)d' % vmd_params
	os.system(vmdstring)


def makeHist(vmd_outfile, outpicklePrefix, normalize = False, delVMDfiles = False):
	try:
		sasa = np.loadtxt(vmd_outfile)
	except ValueError:
		print '%s not found\n' % vmd_outfile
		return
	nframes = len(sasa)
	sasa_min = 0.98 * np.min(sasa)
	sasa_max = 1.02 * np.max(sasa)
	nbins = 50
	delta = (sasa_max - sasa_min)/float(nbins)
	bin_centers = np.zeros([nbins], np.float64)
	bin_vals = np.zeros([nbins], np.float64)
	for i in range(nbins):
		bin_centers[i] = sasa_min + (i + 0.5)*delta
	
	for n in range(nframes):	
		this_sasa = sasa[n,:]
		for this_n in range(len(this_sasa)):
			assignment = int((this_sasa[this_n] - sasa_min)/delta)
			bin_vals[assignment] += 1.0
	
	bin_vals /= nframes
	if normalize:	bin_vals /= (np.sum(bin_vals) * delta)
	if delVMDfiles:	os.remove(vmd_outfile)
	pickle.dump((bin_centers,bin_vals), open(outpicklePrefix + '.pickle', 'w'))
	


def plot_sasa(outpickle, cgtype, ax):
	try:
		(x,y) = pickle.load(open(outpickle, 'r'))
	except IOError:
		print '%s not found\n' % outpickle
	style = {'AA': '-r', 'SP': 'b-', 'SPLD': 'g--', 'LD': 'b:'}
	ax.plot(x,y, style[cgtype], linewidth = 3, label = cgtype)
	ax.legend()



##### Main
data_dir = os.path.abspath('../data/traj')
analysis_dir = os.path.abspath('../data/analysis')
fftypes = ['wca', 'lj']
cgtypes = ['SP', 'SPLD', 'LD', 'AA']
trajtypes = ['CG', 'AA']

## user input
n_mon = int(sys.argv[1]) #chain  length
n_water = int(sys.argv[2]) #number of waters

# filename formats
trajfile_format = 'c%d_%s_%s.lammpstrj'
vmd_outprefix_format = '%s_sasa_%s_%s'
outpicklePrefix_format = '%s_%s_%s_hist1D_SASA_atom'


fig = plt.figure(figsize = (8,4), facecolor = 'w', edgecolor = 'w')
axs = {'LJ': fig.add_subplot(1,2,1), 'WCA': fig.add_subplot(1,2,2)}


for fftype in fftypes:
	for cgtype in cgtypes:
		trajfile = os.path.join(data_dir, trajfile_format % (n_mon, fftype, cgtype))
		if cgtype == 'AA': 
			mon_type = 3
			trajtype = 'AA'
			mon_start = 3*n_water; mon_end = 3*n_water+n_mon - 1
		else:
			mon_type = 0
			trajtype = 'CG'
			mon_start = 0; mon_end = n_mon - 1
		vmd_outprefix = os.path.join(analysis_dir, vmd_outprefix_format % (trajtype, fftype, cgtype))
		outpicklePrefix = os.path.join(analysis_dir, outpicklePrefix_format % (trajtype, fftype, cgtype))
			
		if cgtype == 'AA':
			run_vmd(traj = trajfile, vmd_outprefix = vmd_outprefix, mon_type = mon_type)
			makeHist(vmd_outfile = vmd_outprefix + '.dat', outpicklePrefix = outpicklePrefix)
		plot_sasa(outpickle = outpicklePrefix + '.pickle', cgtype = cgtype, ax = axs[fftype.upper()])


# plot design
for key in axs.keys():
	axs[key].set_xlabel('SASA_atom ' + r'$(\AA^2/atom)$', fontsize = 'large')
	axs[key].set_ylabel('distribution', fontsize = 'large')
	axs[key].xaxis.set_major_locator(MaxNLocator(nbins = 5, prune = 'both'))
	axs[key].yaxis.set_major_locator(MaxNLocator(nbins = 5, prune = 'both'))
if key == 'WCA':
	axs[key].set_title('WCA (superhydrophobic)')

if key == 'LJ':
	axs[key].set_title('LJ (hydrophilic)')

axs['WCA'].set_yticklabels([])
axs['WCA'].set_ylabel('')

plt.subplots_adjust(left = 0.15, bottom = 0.15, wspace = 0)
plt.savefig('sasa.svg', dpi = 300)

plt.show()

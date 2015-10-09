#!/usr/bin/env python

import numpy as np
import os, sys, pickle
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import sasa_mdtraj

# GLOBALS
data_dir = os.path.abspath('../data/traj')
analysis_dir = os.path.abspath('../data/analysis')
pdb_datafiles = {'AA': os.path.join(data_dir, 'c25_AA_data.pdb'), 'CG': os.path.join(data_dir, 'c25_CG_data.pdb')}
vmd_tclscript = os.path.join(os.getcwd(), 'sasa_vmd.tcl')

# Per atom SASA calculation with vmd
def run_vmd(traj, outprefix, mon_start = 0, mon_end = 25, mon_type = 0, trajtype = 'CG'):
	if not os.path.isfile(traj):
		print '%s not found\n' % traj
		return
	
	#trajtype is depreciated; not used
	vmd_params = {'tclscript': vmd_tclscript , 
	      	      'serial_start': mon_start, 'serial_end': mon_end,
	      	      'start': 0, 'end': -1, 'freq': 5,
	      	      'trajfile': traj,
		      'outprefix': outprefix,
		      'mon_type' : mon_type}

	vmdstring = 'vmd -dispdev text -e %(tclscript)s -args %(trajfile)s %(outprefix)s %(serial_start)d %(serial_end)d %(start)d %(end)d %(freq)d %(mon_type)d' % vmd_params
	os.system(vmdstring)


# Per atom SASA calculation with mdtraj
def run_mdtraj(traj, outprefix, mon_start = 0, mon_end = 25, mon_type = 0, trajtype = 'CG'):
	if not os.path.isfile(traj):
		print '%s not found\n' % traj
		return

	datafile = pdb_datafiles[trajtype]
	sasa_mdtraj(trajfile = traj, datafile = datafile, outprefix = outprefix, 
		    serial_start = mon_start, serial_end = mon_end, mon_type = mon_type,
	    	    start = 0, stop = -1, freq = 5)


# Rules for determining params to be passed to SASA calc routines
def getParams(cgtype, method, n_water, n_mon):
	# return in this order (mon_start, mon_end, mon_type, trajtype)
	# add more rules if more software found to calculate sasa
	if not (cgtype == 'AA'): cgtype = 'CG'
	
	ret = {('AA', 'vmd') : (3*n_water+1, 3*n_water+n_mon, 3, 'AA'),
	       ('CG', 'vmd') : (1, n_mon, 0, 'CG'),
	       ('AA', 'mdtraj') : (3*n_water, 3*n_water+n_mon, 3, 'AA'),
	       ('CG', 'mdtraj') : (0, n_mon, 0, 'CG')}

	return ret[(cgtype, method)] 	


# Making histograms of per atom SASA
def makeHist(outfile, outpicklePrefix, normalize = False, delTempfiles = False):
	try:
		sasa = np.loadtxt(outfile)
	except ValueError:
		print 'Weird ValueError'
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
	if delTempfiles:	os.remove(outfile)
	pickle.dump((bin_centers,bin_vals), open(outpicklePrefix + '.pickle', 'w'))
	

# Plotting SASA distribution
def plot_sasa(outpickle, cgtype, ax):
	try:
		(x,y) = pickle.load(open(outpickle, 'r'))
	except IOError:
		print '%s not found\n' % outpickle
		return
	style = {'AA': '-r', 'SP': 'b-', 'SPLD': 'g--', 'LD': 'b:'}
	ax.plot(x,y, style[cgtype], linewidth = 3, label = cgtype)
	ax.legend()



##### Main
fftypes = ['wca', 'lj']
cgtypes = ['SP', 'SPLD', 'LD', 'AA']
trajtypes = ['CG', 'AA']
methods = {'vmd': run_vmd , 
	   'mdtraj': run_mdtraj} #add more functions here later if need be

## user input
n_mon = int(sys.argv[1]) #chain  length
n_water = int(sys.argv[2]) #number of waters

# filename formats
trajfile_format = 'c%d_%s_%s.lammpstrj'
outprefix_base_fmt = '%s_sasa_%s_%s_%s'
outpicklePrefix_base_fmt = '%s_%s_%s_hist1D_SASA_atom_%s'


for method in methods.keys():
	# fig initialization
	fig = plt.figure(figsize = (8,4), facecolor = 'w', edgecolor = 'w')
	axs = {'LJ': fig.add_subplot(1,2,1), 'WCA': fig.add_subplot(1,2,2)}

	for fftype in fftypes:
		for cgtype in cgtypes:
			trajfile = os.path.join(data_dir, trajfile_format % (n_mon, fftype, cgtype))
			mon_start, mon_end, mon_type, trajtype = getParams(cgtype = cgtype, method = method, n_mon = n_mon, n_water = n_water)			
			
			outprefix = os.path.join(analysis_dir, outprefix_base_fmt % (trajtype, fftype, cgtype, method))
			outpicklePrefix = os.path.join(analysis_dir, outpicklePrefix_base_fmt % (trajtype, fftype, cgtype, method))
			
			fnSASA = methods[method]
			fnSASA(traj = trajfile, outprefix = outprefix, mon_start = mon_start, mon_end = mon_end, mon_type = mon_type, trajtype = trajtype)			
			
			makeHist(outfile = outprefix + '.dat', outpicklePrefix = outpicklePrefix, normalize = True)

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
	figname = 'sasa_%s.svg' % method
	plt.savefig(figname, dpi = 300)

plt.show()

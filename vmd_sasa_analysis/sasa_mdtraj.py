#!/usr/bin/env python

def sasa_mdtraj(trajfile, datafile, outprefix, serial_start = 0, serial_end = 25, mon_type = 0,
	    start = 0, stop = -1, freq = 1):

	import numpy as np
	import os
	import mdtraj as md
	
	# note: mdtraj takes sasa radii in nm
	w_rad = 1.4 / 10.
	mon_rad = 2.0934 / 10.
	
	print 'Starting mdtraj for %s...' % outprefix
	traj = md.load_lammpstrj(filename = trajfile, top = datafile, stride = freq, unit_set = 'real')
	
	# selecting polymer atoms based on atom serial number
	poly_indices = range(serial_start, serial_end)
	polymer_sel = traj.atom_slice(atom_indices = poly_indices, inplace = False)
			
	# low-level fix to change atom radius
	# turns out that mdtraj as internal radius conventions
	# by which it automatically sets radius of a P (phosphorus) atom 
	# to 0.16 nm
	# source: https://github.com/mdtraj/mdtraj/blob/master/mdtraj/geometry/sasa.py
	md.geometry.sasa._ATOMIC_RADII['P'] = mon_rad
		
	# run teh sasa calculator
	sasa = md.shrake_rupley(traj = polymer_sel, probe_radius = w_rad, mode = 'atom')
	sasa *= 100 #conversion from nm^2 to A^2

	outfile = outprefix + '.dat'
	np.savetxt(outfile, sasa)


import sys
sys.modules[__name__] = sasa_mdtraj


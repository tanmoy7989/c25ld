#!/usr/bin/env python

import numpy as np
import os, sys

# Globals
curr_dir = sys.argv[1]
c_lens = [15]

# Dir hierarchy
def build_dir(c_len):
	parent_dir = 'aug15_runs_c%d' % c_len
	print 'Building %s ...' % parent_dir
	return os.path.join(curr_dir, parent_dir)
	
# Generate Lammps script for specific chain length
def genLammpsJobs(parent_dir, c_len, n_water= 1700, n_poly = 1):
	print 'Loading stuff into %s ...' % parent_dir
	makeLammpsHome = os.path.expanduser('~/c25ld/makeLammps')
	os.chdir(makeLammpsHome)
	cmdstring = 'python writelammpsjob.py %s %d %d %d 0' % (parent_dir, c_len, n_water, n_poly)
	os.system(cmdstring)
	os.chdir(curr_dir)
	


# Main
for c_len in c_lens:
	p_dir = build_dir(c_len)
	genLammpsJobs(p_dir, c_len)

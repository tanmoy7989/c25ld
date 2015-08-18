#!/usr/bin/env python

import numpy as np
import os, sys

# Globals
trans_dir = sys.argv[1]
curr_month = 'aug15'
c_lens = [15]

# Dir hierarchy
def build_dir(c_len):
	parent_dir = '%s_runs_c%d' % (curr_month, c_len)
	print 'Building %s ...' % parent_dir
	os.mkdir(os.path.join(trans_dir, parent_dir))
	return os.path.join(trans_dir, parent_dir)
	
# Generate Lammps script for specific chain length
def genLammpsJobs(parent_dir, c_len, n_water= 1700, n_poly = 1):
	print 'Loading stuff into %s ...' % parent_dir
	cmdstring = 'python writelammpsjob.py %s %d %d %d 0' % (parent_dir, c_len, n_water, n_poly)
	os.system(cmdstring)
	

# Main
for c_len in c_lens:
	p_dir = build_dir(c_len)
	genLammpsJobs(p_dir, c_len)

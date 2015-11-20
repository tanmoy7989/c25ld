#!/usr/bin/env python

import numpy as np
import os, sys

# Globals
trans_CGdir = sys.argv[1]
trans_AAdir = sys.argv[2]
curr_month = 'aug15'
c_lens = [10, 20, 30, 40, 50, 80]
c25ffdir = os.path.expanduser('~/c25ld/data/cg_ff/feb15_runs_fsw')

# Dir hierarchy
def build_dir(c_len):
	parent_dir = '%s_runs_c%d' % (curr_month, c_len)
	print 'Building %s ...' % parent_dir
	os.mkdir(os.path.join(trans_CGdir, parent_dir))
	CGdir = os.path.join(trans_CGdir, parent_dir)
	AAdir = os.path.join(trans_AAdir, parent_dir)
	return (CGdir, AAdir)
	
# Generate sim script for specific chain length
def genSimJobs(CGdir, AAdir, c_len, n_water= 1700, n_poly = 1):
	print 'Loading stuff into %s ...' % CGdir
	cmdstring = 'python write_MD_jobs.py %d %s %s %s %d %s %s' % (c_len, c25ffdir, CGdir, AAdir, 0, '', '')
	os.system(cmdstring)
	

# Main
for c_len in c_lens:
	CGdir, AAdir = build_dir(c_len)
	genSimJobs(CGdir, AAdir, c_len)

#!/usr/bin/env python

"""
functions that parse out
initial and post-npt box-lengths from
AA trajectories of explicit water simulations
"""

import os, sys
sys.path.append(os.path.expanduser('~/.scripts'))
curr_dir = os.getcwd()

def getInitBoxLen(datafile):
	of = open(datafile, 'r')
	for line in of:
		if line and line.strip().endswith('xhi'):
			l = line.split()
			boxlen = float(l[1]) - float(l[0])
	
	return boxlen
	

def getNPTBoxLen(restartfile, ff_file):
	#if not os.path.isfile(restartfile):
	#	raise IOError('%s not present' % restartfile.split('/')[-1])
	
	
	os.chdir(os.path.dirname(restartfile))
	new_datafile = os.path.join(os.path.dirname(restartfile), 'tmp_res2data.dat')
	cmdstring = 'restart2data.py %s %s %s' % (restartfile, new_datafile, ff_file)
	os.system(cmdstring)
	of = open(new_datafile, 'r')
	for line in of:
		if line.strip().endswith('xhi'):
			l = line.split()
			boxlen = float(l[1]) - float(l[0])
	os.remove(new_datafile)
	return boxlen
	
	
##test
#restartfile = sys.argv[1]
#ff_file = sys.argv[2]
#print getNPTBoxLen(restartfile, ff_file)

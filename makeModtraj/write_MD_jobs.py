#!/usr/bin/env python

import os, sys
import numpy as np
import pickle

#sys.path.append(os.path.expanduser('~/c25ld/analysis'))
from getBoxL import getNPTBoxLen

## Constant values
fftypes = ['lj', 'wca']
cgtypes = ['SP', 'SPLD', 'LD']

sJobIn = file(os.path.expanduser('~/job.template')).read()

## changed this function to get 
## correct BoxLen appropriate to the chain length
## T.S. 08/26/2015
def getBoxL(AAdir, fftype, c_len):
	parent_dir = os.path.join(AAdir, 'unconstrained_%s' % fftype)
	restartfile = os.path.join(parent_dir, 'c%d_unbiased_relax.curr' % c_len)
	ff_file = os.path.join(parent_dir, 'c%d.forcefield' % c_len)

	boxlen = getNPTBoxLen(restartfile, ff_file)
	return boxlen
	

def unconst_job(ffdir, targetdir, AAdir, c_len):
	for fftype in fftypes:
		
		boxlen = getBoxL(AAdir, fftype, c_len)
		
		if not os.path.isdir(os.path.join(targetdir, 'unconstrained_'+fftype)):
			os.mkdir(os.path.join(targetdir, 'unconstrained_' + fftype))
		
		for cgtype in cgtypes:
			d = {}
			this_targetdir = os.path.join(targetdir, 'unconstrained_'+fftype,cgtype)
			if not os.path.isdir(this_targetdir):	os.mkdir(this_targetdir)
			os.system('cp cg_MD.py %s' % this_targetdir)
			
			ffield_file = os.path.join(ffdir, 'c25_%s_%s_sum.txt' % (fftype,cgtype))
			Prefix = 'c%d_%s_%s' % (c_len, fftype, cgtype)
			
			d['JOBNAME'] = 'sim_MD'
			d['CMD'] = 'python cg_MD.py %d %s %s %s %g' % (c_len, ffield_file, Prefix, fftype, boxlen) 
			d['hasMAIL'] = ''
			sJob = os.path.join(this_targetdir, Prefix + '.sh')
			file(sJob, "w").write(sJobIn % d)
			os.system('chmod 777 ' + sJob)
			

def const_job(ffdir, targetdir, biasfiles):
	for fftype in fftypes:
	
		boxlen = getBoxL(ffdir, fftype)
		
		if not os.path.isdir(os.path.join(targetdir, 'constrained_'+fftype)):
			os.mkdir(os.path.join(targetdir, 'constrained_' + fftype))
		
		u_data = np.loadtxt(biasfiles[fftype], skiprows = 0)
		Rg_centers = u_data[:,0]
		forceconst = u_data[:,1]
		
		for cgtype in cgtypes:
			d = {}
			this_targetdir = os.path.join(targetdir, 'constrained_'+fftype,cgtype)
			if not os.path.isdir(this_targetdir):	os.mkdir(this_targetdir)
			os.system('cp cg_umbrella_MD.py %s' % this_targetdir)
			ffield_file = os.path.join(ffdir, 'c25_%s_%s_sum.txt' % (fftype,cgtype))
			
			
			for i in range(len(Rg_centers)):
				Prefix = 'c25_%s_%s_%d' % (fftype, cgtype,i+1)			
				s = 'python cg_umbrella_MD.py %s %s %s %g %g %g' % (ffield_file, Prefix, fftype, boxlen, Rg_centers[i], forceconst[i])
				
				d['JOBNAME'] = 'sim_umbrella'
				d['CMD'] = s
				d['hasMAIL'] = '#'
				sJob = os.path.join(this_targetdir, Prefix + '.sh')
				file(sJob, "w").write(sJobIn % d)
				os.system('chmod 777 ' + sJob)
	


## Input files taken from command prompt
c_len = int(sys.argv[1])
ffdir = os.path.abspath(sys.argv[2])
targetdir = os.path.abspath(sys.argv[3])
AAdir = os.path.abspath(sys.argv[4])
if not AAdir: AAdir = ffdir

#doConstraint = bool(sys.argv[5])
#if doConstraint:
#	biasfile_lj = sys.argv[6]
#	biasfile_wca = sys.argv[7]
#	biasfiles = {'lj': os.path.abspath(biasfile_lj), 'wca': os.path.abspath(biasfile_wca)}
#	const_job(ffdir, targetdir, biasfiles)

unconst_job(ffdir, targetdir, AAdir, c_len)

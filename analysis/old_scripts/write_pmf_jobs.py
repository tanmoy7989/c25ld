#!/usr/bin/env python

import os, sys
job_template = os.path.expanduser('~/job.template')
sJobIn = file(job_template).read()

def get_params(thisdir, trjtype, fftype, N_poly):
	basename = 'c%d' % N_poly
	if trjtype == 'AA':
		trjprefix = basename + 'Rg'
		prefix = 'AA_%s' % fftype
		trjsrc = 'lmp'
	if trjtype == 'CG':
		cgtype = thisdir.split('/')[-1]
		trjprefix = '%s_%s_%s_' % (basename, fftype, cgtype)
		prefix = 'CG_%s_%s' % (fftype, cgtype)
		trjsrc = 'sim'
		
	return {'trjprefix': trjprefix, 'prefix': prefix, 'trjsrc': trjsrc}


## User input
if sys.argv[1]:	AAdir = os.path.abspath(sys.argv[1])
else:	AAdir = ''
if sys.argv[2]:	CGdir = os.path.abspath(sys.argv[2])
else:	CGdir = ''
if sys.argv[3]:	AA_paramfile = os.path.abspath(sys.argv[3])
if sys.argv[4]: CG_paramfile = os.path.abspath(sys.argv[4])
target_dir = os.path.abspath(sys.argv[5])
fftype = sys.argv[6]
if sys.argv[7]:	N_poly = int(sys.argv[7])
else:	N_poly = 25
	
# All-atom jobs
if AAdir:
	dir_fftype = os.path.join(AAdir, 'constrained_%s' % fftype)
	params = get_params(dir_fftype, 'AA', fftype, N_poly)
	sJob = os.path.join(target_dir, 'AA_%s_pmf.sh' % fftype)
	d = {'JOBNAME': 'pmf', 'hasMAIL': '#', 
		 'CMD': 'python pmf_analyse.py %s %s %s %s %s %s' % (dir_fftype, AA_paramfile, 
		 													 params['trjsrc'], params['trjprefix'],
		 													 params['prefix'], target_dir)}
	file(sJob, 'w').write(sJobIn % d)
	os.system('chmod 777 %s' % sJob)

# CG-jobs
if CGdir:
	cgtypes = ['SP', 'SPLD', 'LD']
	for cgtype in cgtypes:
		dir_fftype_cgtype = os.path.join(CGdir, 'constrained_%s' % fftype, cgtype)
		params = get_params(dir_fftype_cgtype, 'CG', fftype, N_poly)
		sJob = os.path.join(target_dir, 'CG_%s_%s_pmf.sh' % (fftype, cgtype))
		d = {'JOBNAME': 'pmf', 'hasMAIL': '#',
			 'CMD': 'python pmf_analyse.py %s %s %s %s %s %s' % (dir_fftype_cgtype, CG_paramfile, 
			 													 params['trjsrc'], params['trjprefix'],
			 													 params['prefix'], target_dir)}
		file(sJob, 'w').write(sJobIn % d)
		os.system('chmod 777 %s' % sJob)

	 												   	
#os.system('cp pmf_analyse.py %s' % target_dir)	 
	 
	 
	 

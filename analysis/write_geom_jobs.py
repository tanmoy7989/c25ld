#/usr/bin/env python

import os, sys
import numpy as np

sJobIn = file('/home/cask0/home/tsanyal/job.template').read()

makeAA = False
makeLD = True

def makeJobDict(fftype, AAdir, CGdir, target_dir):
	
	T_AA = {'name' : os.path.join(AAdir, 'unconstrained_%s' % fftype, 'c25_unbiased.lammpstrj.gz'),
		'type' : 'AA_%s' % fftype,
	  	'src' : 'lmp',
	  	'target_dir': target_dir}

	T_SP = {'name' : os.path.join(CGdir, 'unconstrained_%s' % fftype, 'SP', 'c25_%s_SP.lammpstrj.gz' % fftype),
	  'type' : 'CG_%s_SP' % fftype,
	  'src' : 'sim',
	  'target_dir': target_dir}
	  
	
	T_SPLD = {'name' : os.path.join(CGdir, 'unconstrained_%s' % fftype, 'SPLD', 'c25_%s_SPLD.lammpstrj.gz' % fftype),
	  'type' : 'CG_%s_SPLD' % fftype,
	  'src' : 'sim',
	  'target_dir': target_dir}
	  
	if makeLD:
		T_LD = {'name' : os.path.join(CGdir, 'unconstrained_%s' % fftype, 'LD', 'c25_%s_LD.lammpstrj.gz' % fftype),
	  	'type' : 'CG_%s_LD' % fftype,
	  	'src' : 'sim',
	  	'target_dir': target_dir}
	else:
		T_LD = {}
	  
	  
	return [T_AA, T_SP, T_SPLD, T_LD]
	  


def makeInputFile(JobDict):
	
	T = JobDict
	InputFile = {}
	
	T_AA = T[0]
	T_SP = T[1]
	T_SPLD = T[2]
	T_LD = T[3]
	
	InputFile['AA'] = ('AA_%s_input.txt' % fftype, T_AA)
	InputFile['CG_SP'] = ('CG_%s_SP_input.txt' %fftype, T_SP)
	InputFile['CG_SPLD'] = ('CG_%s_SPLD_input.txt' %fftype, T_SPLD)
	if T_LD:
		InputFile['CG_LD'] = ('CG_%s_LD_input.txt' %fftype, T_LD)

	for key in InputFile.keys():
		v = InputFile[key]
		filename = os.path.join(target_dir, v[0])
		print filename
		file(filename, 'w').write(str(v[1]))
		
	return InputFile
	
	


def makeJobs(InputFiles, sJobIn, target_dir, fftype):
	d = {'JOBNAME' : 'analysis', 'CMD': ''}
	for key in InputFiles.keys():
		filename = os.path.join(target_dir, InputFiles[key][0])
		sJob = filename.split('/')[-1].split('_input.txt')[0]
		sJob = os.path.join(target_dir, sJob + '_analysis.sh')
		d['CMD'] = 'python geom_analyse.py %s' % filename
		file(sJob, 'w').write(sJobIn % d)
		os.system('chmod 777 %s' % sJob)



####### MAIN #######

AAdir = os.path.abspath(sys.argv[1])
CGdir = os.path.abspath(sys.argv[2])
target_dir = os.path.abspath(sys.argv[3])

os.system('cp geom_analyse.py %s ' % target_dir)

fftypes = ['lj', 'wca']
for fftype in fftypes:

	jobDict = makeJobDict(fftype, AAdir, CGdir, target_dir)
	InputFiles = makeInputFile(jobDict)
	makeJobs(InputFiles, sJobIn, target_dir, fftype)
		

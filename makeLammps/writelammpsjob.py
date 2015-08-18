#!/usr/bin/env python

import os, sys
import numpy as np
import time

sLammpsIn = file("./templates/lammpsin.template").read()
sColvarsIn = file("./templates/colvar.in.template").read()
sJobIn 	  = {'constrained': file("./templates/const_job.template").read(),
	     'unconstrained' : file("./templates/unconst_job.template").read()}

autoSubmit = False
curr_dir = os.getcwd()

## Umbrella sampling jobs

def constrained_job(target_dir, biasdata, datafileprefix):
	for i in range(len(biasdata)):
		d_in = {"DATAFILEPREFIX": datafileprefix, "RUNPREFIX" :datafileprefix+'Rg'+ str(i), "FIX_STATE" : '',
			"RELAXSTEPS":500000, "EQUILSTEPS": 1000000, "PRODSTEPS": 8000000, "STEPFREQ": 100, "RESTARTFREQ": 1000}
		d_colvars = {"LBOUND": 3.0, "UBOUND": 10.0, "Rg0" : biasdata[i,0], "FConst": biasdata[i,1], 'STEPFREQ': 100}
		d_job = {"RUNPREFIX" : d_in["RUNPREFIX"],
				"NMON": int(datafileprefix.split('c')[1])}
		
		sLammps =  os.path.join(target_dir, d_in["RUNPREFIX"] + ".in")
		sColvars = os.path.join(target_dir, d_in["RUNPREFIX"] + ".colvar.input")
		sJob = os.path.join(target_dir, d_in["RUNPREFIX"] + ".sh")
		
		file(sLammps, "w").write(sLammpsIn % d_in)
		file(sColvars,"w").write(sColvarsIn % d_colvars)
		file(sJob, "w").write(sJobIn['constrained'] % d_job)
		os.system("chmod +x " + sJob)
		if autoSubmit:
			os.chdir(target_dir)
			os.system('smartsubmit ' + sJob)
			os.chdir(curr_dir)


## Unbiased MD jobs	
def unconstrained_job(target_dir, datafileprefix):
	d_in = {"DATAFILEPREFIX": datafileprefix, "RUNPREFIX" : datafileprefix+'_unbiased' , "FIX_STATE" : '#', 
		"RELAXSTEPS": 500000, "EQUILSTEPS": 1000000, "PRODSTEPS": 26000000, "STEPFREQ": 100, "RESTARTFREQ": 1000}
	d_job = {"RUNPREFIX" : d_in["RUNPREFIX"],
			 "NMON": int(datafileprefix.split('c')[1])}
		
	sLammps =  os.path.join(target_dir, d_in["RUNPREFIX"] + ".in")
	sJob = os.path.join(target_dir, d_in["RUNPREFIX"] + ".sh")
		
	file(sLammps, "w").write(sLammpsIn % d_in)
	file(sJob, "w").write(sJobIn['unconstrained'] % d_job)
	os.system("chmod +x " + sJob)
	if autoSubmit:
		os.system(target_dir)
		os.system('smartsubmit ' + sJob)
		os.chdir(curr_dir)
	



### MAIN JOB SCEHDULER SCRIPT ######

## User Input
trajdir = sys.argv[1]
n_mon = int(sys.argv[2])
n_water = int(sys.argv[3])
n_poly = int(sys.argv[4])
doConstrain = bool(int(sys.argv[5]))
if doConstrain:
	biasfile_lj = sys.argv[6]
	biasfile_wca = sys.argv[7]

datafileprefix = 'c%d' % n_mon # changed according to chain length of polymer

if not os.path.isfile(os.path.join(trajdir, datafileprefix+'.data')):
	print "Creating data file..."
	time.sleep(2)
	cmdstring = 'python genLammpsData.py 298.0 %d %d %d %s %s' % (n_water, n_mon, n_poly, datafileprefix, trajdir)
	os.system(cmdstring)
	
print "Creating directory structures..."
if doConstrain:
	dirs = ['constrained_lj', 'constrained_wca', 'unconstrained_lj', 'unconstrained_wca']
else:
	dirs = ['unconstrained_lj', 'unconstrained_wca']
target_dirs = []
[target_dirs.append(os.path.join(trajdir, this)) for this in dirs]

print "Writing jobs..."
for this_dir in target_dirs:
	os.mkdir(this_dir)
	os.system('cp ' + os.path.join(trajdir, datafileprefix+'.data') + ' ' +  this_dir)
	dirname = dirs[target_dirs.index(this_dir)]
	
	if dirname=='constrained_lj':
		biasdata_lj = np.loadtxt(biasfile_lj)
		os.system('cp ./templates/lj_forcefield.template ' + os.path.join(this_dir, datafileprefix+'.forcefield'))
		constrained_job(this_dir, biasdata_lj,datafileprefix)
		
	if dirname=='constrained_wca':
		biasdata_wca = np.loadtxt(biasfile_wca)
		os.system('cp ./templates/wca_forcefield.template ' + os.path.join(this_dir, datafileprefix+'.forcefield'))
		constrained_job(this_dir, biasdata_wca, datafileprefix)
		
	if dirname=='unconstrained_lj':
		os.system('cp ./templates/lj_forcefield.template ' + os.path.join(this_dir, datafileprefix+'.forcefield'))
		unconstrained_job(this_dir, datafileprefix)
		
	if dirname=='unconstrained_wca':
		os.system('cp ./templates/wca_forcefield.template ' + os.path.join(this_dir, datafileprefix+'.forcefield'))
		unconstrained_job(this_dir, datafileprefix)
		

	

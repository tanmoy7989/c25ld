#/usr/bin/env python

"""	
Contains parsers for extracting data from potential and colvar logs for 
sim and lammps logs.
"""

import numpy as np
import os

Trajdir = ''
Prefix = ''
nstages = ''
trjsrc = 'lmp'

# Umbrella log parser
def parse_Ulog():
	print "Parsing colvar traj logs"
	E_u = []
	Rg =  []
	for stage in range(nstages):
		if trjsrc == 'lmp':
			st = stage
			## determine the start of prod step counter in colvars-traj file
			## added 05/26/2015 by TS
			logfilename = os.path.join(Trajdir, Prefix + str(st)+ ".log")
			lines = file(logfilename, 'r').readlines()
			start = [lines.index(line) for line in lines if line.startswith('Step Temp Atoms PotEng polyPE_t Rg')][0]
			start += 2
			prod_step_start = int(lines[start].split()[0])

		elif trjsrc =='sim':
			st = stage + 1
			prod_step_start = 100
		
		print st
		filename = os.path.join(Trajdir, Prefix + str(st) +".colvars.traj")
		rdata = np.loadtxt(filename, skiprows = 1, comments = '#')
		colvar_step_start = np.where(rdata == prod_step_start)[0][0]
		Rg.append(rdata[colvar_step_start:,1])
		E_u.append(rdata[colvar_step_start:,2])
	
	Rg = np.array(Rg)
	E_u = np.array(E_u)
	return (E_u, Rg)


# Potential log parsers
def parse_PElog():
	tokens = {'start': 'Step Temp Atoms PotEng polyPE_t Rg',
			  'stop': 'Loop time of', 
			  'break' : 'colvars'}
	print "Parsing potential energy logs..."
	E_pe = []
	for stage in range(nstages):
		if trjsrc == 'lmp':
			st = stage
		elif trjsrc =='sim':
			st = stage + 1
		filename = os.path.join(Trajdir, Prefix + str(st) + ".log")
		
		if trjsrc == 'lmp':
			lines = open(filename, 'r').readlines()
			start = [lines.index(line) for line in lines if line.startswith(tokens['start'])][0] + 1
			stop = [lines.index(line) for line in lines if line.startswith(tokens['stop'])][-1] 
			rdata = []
			for index, line in enumerate(lines[start:stop]):
				if line.startswith(tokens['break']):	continue
				rdata.append(float(line.split()[-2].strip()))
		elif trjsrc == 'sim':
			rdata = np.loadtxt(filename, skiprows = 1)[:,2]
		E_pe.append(rdata)
	
	E_pe = np.array(E_pe)
	return E_pe

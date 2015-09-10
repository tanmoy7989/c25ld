#!/usr/bin/env python


import numpy as np
import os, sys, pickle
import sim

## Assumption: required chem_data and makesys are in the 
## same directory
import makesys
from chem_data import getModelData

d = getModelData()
poly_dict = d['poly_dict']
TempSet = d['TempSet']

NStep = {'Stepsize': 0.003, 'Equil': 1000000, 'Prod': 10000000, 'Freq': 100}


def parseTraj(lammpstraj):
	"""
	1)Pickles trajectory to save time in reading frames
	2) Parses out system name from trajectory name
	3) Parses out length of box
	4) Returns ids of atoms that are to be mapped in the CG model
	"""
	
	monomer_atomtype = '3'
	if not os.path.isfile(lammpstraj + ".pickle"):
		print "Pickling trajectory...."
		Trj = sim.traj.lammps.Lammps(lammpstraj)
		pickle.dump(Trj, open(lammpstraj + '.pickle', 'w'))
	else:
		print 'Loading pre-Pickled trajectory...'
		Trj = pickle.load(open(lammpstraj + '.pickle', 'r'))
	Pos = Trj[0]
	boxlen = Trj.FrameData["BoxL"]
	mappedatoms = []
	for (i, AtomName) in enumerate(Trj.AtomNames):
		if AtomName == monomer_atomtype:
			mappedatoms.append(i)
	if not len(mappedatoms) == poly_dict["N_mon"] * poly_dict["N_poly"]:
		raise ValueError("Number of atoms to map in trajectory does not match that in polydict")
	
	return (Trj, boxlen, mappedatoms)



## Parse input reference trajectory
lammpstraj = sys.argv[1]
fftype = sys.argv[2]
BasePrefix = sys.argv[3]
DataLogPrefix = BasePrefix

(Trj, boxlen, mappedatoms)  = parseTraj(lammpstraj)
of = open(DataLogPrefix + '.dat', 'w')
of.write('#Cutoff\t Srel\n')
of.close()

## Print details
print "Reference Lammps Trajectory =  ", lammpstraj.split('/')[-1]
print "System Temperature  = ", TempSet
print "Number of polymer molecules = ", poly_dict["N_poly"]
print "Number of monomers = ", poly_dict["N_mon"]
print "Length of simulation box = ", boxlen

## Start cutoff loop
cutoff_list = np.arange(3.0, 8.0,0.1)
for i, cutoff in enumerate(cutoff_list):
	
	## Misc
	print 'Local Density Cutoff = %d' % cutoff
	Prefix = BasePrefix + '_%d' % i
	
	## Create the energy-minimized System object
	Sys = makesys.MakeSys(fftype = fftype, BoxL = boxlen, Prefix = Prefix, LDCut = cutoff)
	Int = Sys.Int

	#shortcuts to potentials
	Plj = Sys.ForceField[0]
	Pspline = Sys.ForceField[1]
	Pld = Sys.ForceField[2]

	#change to MD
	Int.Method = Int.Methods.VVIntegrate
	Int.Method.TimeStep = NStep['Stepsize']
	Int.Method.Thermostat = Int.Method.ThermostatLangevin
	Int.Method.LangevinGamma = 0.01
	Sys.TempSet = TempSet

	## Mapping function
	Map = sim.atommap.PosMap()
	for (i, a) in enumerate(Sys.Atom):
	    Map += [sim.atommap.AtomMap(Atoms1 = mappedatoms[i], Atom2 = a)]

	## Spline treatment
	Pspline.EneSlopeInner = None

	## Set up Srel optimizer
	Sys.TempSet = TempSet
	Opt = sim.srel.OptimizeTrajClass(Sys, Map, Traj = Trj, SaveLoadArgData = True, FilePrefix = Prefix)
	Opt.FilePrefix = Prefix 
	Opt.MinReweightFrac = 0.20
	sim.srel.optimizetraj.PlotFmt = 'svg'
	
	## Lammps settings
	sim.export.lammps.LammpsExec = '/home/cask0/home/tsanyal/tanmoy_lammps/lammps-15May15/src/lmp_ZIN'
	sim.export.lammps.InnerCutoff = 0.01

	## Start srel optimization
	Pspline.SetParam(Knots = 0.)
	Pld.SetParam(Knots = 0.)
	Sys.ForceField.Update()
	Opt.RunConjugateGradient(StepsEquil = NStep['Equil'], StepsProd = NStep['Prod'], StepsStride = NStep['Freq']) 
	
	## Log data
	lines = open(Prefix + '_log.txt', 'r').readlines()
	srel = float(lines[-1].split()[2]) 
	of = open(DataLog + '.dat', 'a')
	of.write('%g\t %g\n' % (cutoff, srel))
	of.close()

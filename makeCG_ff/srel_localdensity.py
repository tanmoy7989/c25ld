#!/usr/bin/env python


import numpy as np
import os, sys, pickle
import sim

sys.path.append(os.path.expanduser('~/c25ld'))

import makesys
from chem_data import getModelData

d = getModelData()
poly_dict = d['poly_dict']
cut_dict = d['cut_dict']
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
Prefix = sys.argv[3]

(Trj, boxlen, mappedatoms)  = parseTraj(lammpstraj)
print "Reference Lammps Trajectory =  ", lammpstraj.split('/')[-1]
print "System Temperature  = ", TempSet
print "Number of polymer molecules = ", poly_dict["N_poly"]
print "Number of monomers = ", poly_dict["N_mon"]
print "Length of simulation box = ", boxlen
	

## Create the energy-minimized System object
Sys = makesys.MakeSys(fftype = fftype, BoxL = boxlen, Prefix = Prefix)
Int = Sys.Int

#shortcuts to potentials
Plj = Sys.ForceField[0]
Pspline = Sys.ForceField[1]
Pbond = Sys.ForceField[2]
Pangle = Sys.ForceField[3]
Pld = Sys.ForceField[4]


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
sim.srel.optimizetraj.PlotFmt = 'svg'

## Output initial histograms
#Opt.MakeModTraj(StepsEquil=10000, StepsProd=1000000, StepsStride=100)
#Opt.OutputModHistFile()
#Opt.OutputPlot()

## Lock all potentials in place
def ResetPotentials():
	print "Resetting all potentials to base case and freezing."
	Plj.SetParam(Epsilon = poly_dict["LJEpsilon"], Sigma = poly_dict["LJSigma"])
	Pspline.SetParam(Knots = 0.)
	Pld.SetParam(Knots = 0.)
	for ptype in Sys.ForceField:
		ptype.FreezeParam()


## Srel minimization
opt_cases = ["SP", "SPLD", "LD"]

for i, case in enumerate(opt_cases):
    Opt.Reset()
    Opt.FilePrefix = Prefix + "_" + case
    print "\n\nOptimizing for the case: ", case
    
    if case == "SP":
        ResetPotentials()
        Plj.SetParam(Epsilon = 0., Sigma = 1.e-10)
        Pspline.UnfreezeParam()
        print "Turned off LJ potential."
        print "Unfrozen potential = SP"

    elif case == "SPLD":
        Pld.UnfreezeParam()
        print "Unfrozen potentials = SP, LD"

    elif case == "LD":
        ResetPotentials()
        Pld.UnfreezeParam()
        print "Unfrozen potential = LD"
    	
       
    Sys.ForceField.Update()
    Opt.RunConjugateGradient(StepsEquil = NStep['Equil'],
                             StepsProd = NStep['Prod'],
                             StepsStride = NStep['Freq']) 

print "Srel minimization for different cases finished"

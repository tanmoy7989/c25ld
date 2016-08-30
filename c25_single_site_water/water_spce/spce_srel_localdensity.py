#!/usr/bin/env python

import numpy as np
import os, sys, pickle
import sim

np.random.seed(12345)

SPKnots = [1.0178e+02 , 9.1777e+01 , 8.1777e+01 , 7.1777e+01 , 6.1777e+01 , 5.1777e+01 , 4.1777e+01 , 3.1777e+01 , 2.1777e+01 , 1.1777e+01 , 1.7770e+00 , -3.2263e-01, 4.8323e-01 , 9.2257e-01 , 7.9679e-01 , 5.5024e-01 , 3.6377e-01 , 1.7821e-01 , 1.1474e-02 , -7.1876e-02, -4.8802e-02, 4.8243e-02 , 1.3728e-01 , 8.9832e-02 , 4.6524e-02 , 2.8392e-02 , 1.4874e-02 , 2.2107e-02 , 9.5437e-03 , 3.4601e-02 ]

LDKnots =  [2.3252e-01 , 2.3252e-01 , 2.3252e-01 , 2.3252e-01 , 2.3252e-01 , 2.3252e-01 , 2.3006e-01 , 2.2260e-01 , 2.0294e-01 , 1.1418e-01 , 0.0000e+00 , -1.8159e-01, -3.9348e-01, -6.2700e-01, -8.2742e-01, -8.2742e-01, -8.2742e-01, -8.2742e-01, -8.2742e-01, -8.2742e-01]

TempSet = 298.0
cut_dict = {"SPCut": 7.5, "LDCut": 3.5, "LDLowerCut" : 3.0}            
N_mol = 500

def parseTraj(lammpstraj):
	monomer_atomtype = '1'
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
	if not len(mappedatoms) == N_mol:
		raise ValueError("Number of atoms to map in trajectory does not match that in spce_dict")
	
	return (Trj, boxlen, mappedatoms)


## Parse input reference trajectory
lammpstraj = sys.argv[1]
(Trj, boxlen, mappedatoms)  = parseTraj(lammpstraj)
print "Reference Lammps Trajectory =  ", lammpstraj.split('/')[-1]
print "System Temperature  = ", TempSet
print "Number of water molecules = ", N_mol
print "Length of simulation box = ", boxlen
	
## Making the world and the system
sysname = 'spce'
atomtype = sim.chem.AtomType("W", Mass = 18.0142, Charge = 0.0, Color = (1,0,0), Radius = 0.5)
moltype = sim.chem.MolType("W", [atomtype])
world = sim.chem.World([moltype], Dim = 3, Units = sim.units.AtomicUnits)
Sys = sim.system.System(world, Name = sysname)
for i in range(N_mol): Sys += moltype.New()

# Assigning box length
Sys.BoxL[:] = boxlen

## Desiging the potentials
Pspline = sim.potential.PairSpline(Sys, Filter = sim.atomselect.NonbondPairs,
                                   Cut = cut_dict["SPCut"], NKnot = len(SPKnots), Label = "SP")


LDFilter = sim.atomselect.PolyFilter([atomtype, atomtype], Ordered = True)
LDRhoMin = 0.0 ; LDRhoMax = 20.0
Plocaldensity = sim.potential.LocalDensity(Sys, Filter = LDFilter,
                                            Cut = cut_dict["LDCut"], LowerCut = cut_dict["LDLowerCut"],
                                            RhoMin = LDRhoMin, RhoMax = LDRhoMax, NKnot = len(LDKnots),
                                            Label = "LD")
                                                                                    
## Spline treatment
Sys.ForceField.SetSplineTreatment( NonbondEneSlope = 40., BondEneSlope = 10., AngleEneSlope = 50.)
Sys.ForceField.extend([Pspline, Plocaldensity])


## Integrator settings
Int = Sys.Int
Int.Method = Sys.Int.Methods.VVIntegrate
for P in Sys.ForceField:    
	P.Arg.SetupHist(NBin = 10000, ReportNBin = 100)
Sys.Load()

## Initialization
sim.system.init.positions.CubicLatticeFill(Sys, Random = 0.1)
sim.system.init.velocities.Canonical(Sys, Temp = TempSet)

#run for a bit to minimize and remove clashes
#Int.Method = Int.Methods.VVQuench
#Int.Run(1000, "Minimizing")

#change to MD
Int.Method = Int.Methods.VVIntegrate
Int.Method.Thermostat = Int.Method.ThermostatLangevin
Int.Method.LangevinGamma = 0.02
Sys.TempSet = TempSet

## Mapping function
Map = sim.atommap.PosMap()
for (i, a) in enumerate(Sys.Atom):
    Map += [sim.atommap.AtomMap(Atoms1 = mappedatoms[i], Atom2 = a)]
    
## Set up Srel optimizer
Sys.TempSet = TempSet
Opt = sim.srel.OptimizeTrajLammpsClass(Sys, Map, Traj = Trj, SaveLoadArgData = True, FilePrefix = sysname)
sim.srel.optimizetraj.PlotFmt = 'svg'

## Lammps settings
sim.export.lammps.LammpsExec = '/home/cask0/home/tsanyal/software/tanmoy_lammps/lammps-15May15/src/lmp_ZIN'

## Freeze all potentials
print 'All potentials frozen'
Pspline.SetParam(Knots = 0.) ; Plocaldensity.SetParam(Knots = 0.)
Pspline.FreezeParam()
Plocaldensity.FreezeParam()

print Sys.ForceField.ParamString
			
## Srel minimization
cases = ['SP', 'SPLD']
for case in cases:
	print "Srel optimization for : ", case
	Prefix = sysname + "_" + case
	Opt.Reset()
	Opt.FilePrefix = Prefix

	if case == 'SP':
	   Pspline.UnfreezeParam()
	   Pspline.SetParam(Knots = SPKnots)
	   print 'Unfreezing Spline potential'
		
	if case == 'SPLD':
	   Plocaldensity.UnfreezeParam()
	   Plocaldensity.SetParam(Knots = LDKnots)
	   print 'Unfreezing Spline and Local Density potentials'
		
	Sys.ForceField.Update()
	Opt.RunConjugateGradient(StepsEquil = 500000,
                             StepsProd = 1000000,
                             StepsStride = 100)


print "Srel minmization for different cases finished"

#!/usr/bin/env python

import numpy as np
import os, sys, pickle
import sim

## dependencies
sys.path.append('/home/cask0/home/tsanyal/srel_rc')
import alg
import pickleTraj
import parse_potential as pp

NStep_methane = {'Stepsize': 0.003, 'Min': 1000, 'Equil': 2000000, 'Prod': 20000000, 'Freq': 100}
NStep_polymer = {'Stepsize': 0.003, 'Min': 1000, 'Equil': 2000000, 'Prod': 20000000, 'Freq': 100}

def parseTraj(lammpstraj):
	monomer_atomname = '3'
	Trj = pickleTraj(lammpstraj)
	boxlen = Trj.FrameData["BoxL"]
	mappedatoms = []
	for (i, AtomName) in enumerate(Trj.AtomNames):
		if AtomName == monomer_atomname:
			mappedatoms.append(i)
	return (Trj, boxlen, mappedatoms)


#user input
fftype = sys.argv[1]
isPolymer = bool(int(sys.argv[2]))
hasLJ = bool(int(sys.argv[3]))
N_mon = int(sys.argv[4])
N_water = int(sys.argv[5])
LammpsTraj = sys.argv[6] 
Prefix = sys.argv[7]

LDCut = float(sys.argv[8])
if len(sys.argv) > 9: checkpt = sys.argv[9]
else: checkpt = 'start'

#MD iteration settings
if isPolymer:
    NStep = NStep_polymer
else:
    NStep = NStep_methane
    
#system parameters
SysData = alg.makeSysData(isPolymer = isPolymer, hasLJ = hasLJ, N_mon = N_mon, N_water = N_water)
TempSet = 298.0
Delta = 1.2

#AA traj parameters
Trj, BoxL, mappedatoms = parseTraj(LammpsTraj)

#create the system object
Sys = alg.makeSys(LDCut = LDCut, SysData = SysData, fftype = fftype, BoxL = BoxL, Prefix = Prefix, paramstring = None, Delta = 1.2)

#initialize and energy minimize
sim.system.init.positions.CubicLattice(Sys, Random = 0.1)
sim.system.init.velocities.Canonical(Sys, Temp = TempSet)
Int = Sys.Int
Int.Method = Int.Methods.VVQuench
Int.Run(NStep['Min'])

#change to MD
Int.Method = Int.Methods.VVIntegrate
Int.Method.Thermostat = Int.Method.ThermostatLangevin
Int.Method.LangevinGamma = 0.01
Sys.TempSet = TempSet

#mapping function
Map = sim.atommap.PosMap()
for (i, a) in enumerate(Sys.Atom):
    Map += [sim.atommap.AtomMap(Atoms1 = mappedatoms[i], Atom2 = a)]

#extract the potentials into variables preserving correct sequence
Plj = None; Pspline = None; Pbond = None; Pangle = None; Plocaldensity = None
if hasLJ:
    Plj = Sys.ForceField[0]
    Pspline = Sys.ForceField[1]
else:
    Pspline = Sys.ForceField[0]

if isPolymer:
    Pbond = Sys.ForceField[-3]
    Pangle = Sys.ForceField[-2]

Plocaldensity = Sys.ForceField[-1]

#spline treatment
if not Pspline is None: 
    if fftype == 'lj': Sys.ForceField.SetSplineTreatment(NonbondEneSlope = 50., BondEneSlope = 10., AngleEneSlope = 60.)
    else: Pspline.EneSlopeInner = None

#set up Srel optimizer
Sys.TempSet = TempSet
Opt = sim.srel.OptimizeTrajClass(Sys, Map, Traj = Trj, SaveLoadArgData = True, FilePrefix = Prefix)
sim.srel.optimizetraj.PlotFmt = 'svg'
Opt.MinReweightFrac = 0.2

#lock all potentials in place
def ResetPotentials():
	print "Resetting all potentials to base case and freezing."
	if not Pspline is None: Pspline.SetParam(Knots = 0.)
	if not Plocaldensity is None:  Plocaldensity.SetParam(Knots = 0.)
	if not Plj is None: Plj.SetParam(Epsilon = 0.)
	for ptype in Sys.ForceField:
		ptype.FreezeParam()

#restarting a failed cg model
def canRestartFrom(checkpt = 'SP'):
    return os.path.isfile('%s_%s_sum.txt' % (Prefix, checkpt))
    
def Restart(checkpt = 'SP'):
    if checkpt == 'start':
        return
    else:
        print 'Restarting from after ', checkpt
        checkptFile = Prefix + '_%s_sum.txt' % checkpt
        paramstring = pp.parseParamString(checkptFile)
        Sys.ForceField.SetParamString(paramstring)
        Sys.ForceField.Update()
        for p in Sys.ForceField: p.FreezeParam()
        
        if checkpt == 'SP': 
            Pspline.UnfreezeParam()
        
        if checkpt == 'SPLD':
            Pspline.SetParam(Knots = 0.)
            Plocaldensity.UnfreezeParam()
            if not Plj is None: Plj.SetParam(Epsilon = SysData['poly_dict']['LJEpsilon'])
        
        if checkpt == 'LD':
            Plocaldensity.UnfreezeParam()
            

print Sys.ForceField.ParamString()


#srel minimization
opt_cases = ["SP", "SPLD", "LD"]
if opt_cases.__contains__(checkpt):
    x = opt_cases.index(checkpt)
    opt_cases = opt_cases[x:]
for i, case in enumerate(opt_cases):
    Opt.Reset()
    Opt.FilePrefix = Prefix + '_' + case
    print "\n\nOptimizing for the case: ", case
    if case == "SP":
        if canRestartFrom(case):
            Restart(case)
            continue
        else:
            ResetPotentials()
            Pspline.UnfreezeParam()
            print "Unfreezing SP..."

    if case == "SPLD":
        if canRestartFrom(case):
            Restart(case)
            continue
        else:
            Plocaldensity.UnfreezeParam()
            print "Unfreezing LD..."
    
    if case == "LD":
        if canRestartFrom(case):
            Restart(case)
        else:
            if not Plj is None: Plj.SetParam(Epsilon = SysData['poly_dict']['LJEpsilon'])
            Plocaldensity.UnfreezeParam()
            print "Unfreezing LD..."
        
        
    Sys.ForceField.Update()
    print Sys.ForceField.ParamString()
    Opt.RunConjugateGradient(StepsEquil = NStep['Equil'], StepsProd = NStep['Prod'], StepsStride = NStep['Freq']) 

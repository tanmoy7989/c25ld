#!/usr/bin/env python

import numpy as np
import os, sys
import sim

## dependencies
sys.path.append('/home/cask0/home/tsanyal/srel_rc')
import alg
import pickleTraj
import parse_potential as pp

NStep_methane = {'Stepsize': 0.003, 'Min': 1000, 'Equil': 2000000, 'Prod': 30000000, 'Freq': 100}
NStep_polymer = {'Stepsize': 0.003, 'Min': 1000, 'Equil': 2000000, 'Prod': 30000000, 'Freq': 100}

#user input
fftype = sys.argv[1]
isPolymer = bool(int(sys.argv[2]))
hasLJ = bool(int(sys.argv[3]))
N_mon = int(sys.argv[4])
N_water = int(sys.argv[5])
LammpsTraj = sys.argv[6] 
Prefix = sys.argv[7]
cgffPrefix = sys.argv[8]
LDCut = float(sys.argv[9])

#MD iteration settings
if isPolymer:
    NStep = NStep_polymer
else:
    NStep = NStep_methane

#system parameters
SysData = alg.makeSysData(isPolymer = isPolymer, hasLJ = hasLJ, N_mon = N_mon, N_water = N_water)
TempSet = 298.0
Delta = 1.2

#CG forcefield parameters
cg_ff_file = cgffPrefix + '_sum.txt'
paramstring = pp.parseParamString(cg_ff_file)

#boxlength extracted from AA traj
BoxL = pickleTraj(LammpsTraj).FrameData['BoxL']

#create the system object
Sys = alg.makeSys(LDCut = LDCut, SysData = SysData, fftype = fftype, BoxL = BoxL, Prefix = Prefix, paramstring = paramstring, Delta = 1.2)

#initialize system
sim.system.init.positions.CubicLattice(Sys, Random = 0.1)
sim.system.init.velocities.Canonical(Sys, Temp = TempSet)

##different MD steps
Int = Sys.Int

#energy minimization
print 'Minimizing energy...'
Int.Method = Int.Methods.VVQuench
Int.Run(NStep['Min'])
				
#equilibration
print "Equilibrating..."
Int.Method = Int.Methods.VVIntegrate
Int.Method.TimeStep = NStep['Stepsize']
Sys.Measures.Reset()
Sys.Measures.VerboseOutput(StepFreq = NStep["Freq"])
Sys.TempSet = TempSet
Int.Reset()
Int.Run(NStep['Equil'], 'Equilibrating')

#production
print "Production..."
Int.Method.Thermostat = Int.Method.ThermostatNoseHoover
Sys.Measures.Reset()
#add a Lammps format trajectory dump action
Trj = sim.traj.lammps.LammpsWrite(FileName = Prefix + '.lammpstrj.gz', MinImage = True)
Trj.AddAction(Int, StepFreq = NStep['Freq'])
Int.Reset()  
Int.Run(NStep['Prod']) 

#cleanup
Trj.DelAction(Int)
Trj.Close()

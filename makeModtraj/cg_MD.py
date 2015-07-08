
#!/usr/bin/env python

import numpy as np
import os
import sim
import sys

useVis = False # Flag for visual-python tool

sys.path.append(os.path.expanduser('~/c25ld'))
import makesys
from chem_data import getModelData

d = getModelData()
poly_dict = d['poly_dict']
cut_dict = d['cut_dict']
TempSet = d['TempSet']

## Iteration numbers
NStep = {"Stepsize": 0.003, "StepFreq": 1000, "Equil": 1000000, "Prod": 20000000}
 


def parseFF(filename):
    """
    parses out forcefield parameters from given filename
    """
    with open(filename, "r") as of: lines = of.readlines()
    of.close()
    start = [lines.index(line) for line in lines if line.startswith(">>>")][0]
    stop =  [lines.index(line) for line in lines if line.find("POTENTIAL SUMMARIES")>0][0] - 1
    s = ''
    for line in lines[start:stop]:
	s+= str(line) 
    return s


def writeToFile(PElog, Sys):
    """ 
    Integrator action that logs total potential
    and umbrella run details
    """		
    def PrintHead(Sys):
        Sys.Measures.Reset()
        with open(PElog, "w") as fPE:
        	fPE.write("%7s " % "Step" + Sys.Measures.ReportHead()+"\n")

    def PrintMeasures(Sys):
        with open(PElog, "a") as fPE:
        	fPE.write("%7d " % Sys.Int.StepNum + Sys.Measures.ReportVals() + "\n")


    action = sim.integrate.Action(StepFreq = NStep["StepFreq"], Fn = PrintMeasures, InitFn = PrintHead)
    return action


## Input from user
## Input format:  python cg_MD.py ffield_file Prefix fftype 
ffield_file = sys.argv[1]
Prefix = sys.argv[2]
fftype = sys.argv[3]
boxlen = sys.argv[4]
cgtype = Prefix.split('_')[-1]

## Generating an energy minimized system
Sys = makesys.MakeSys(fftype, BoxL = boxlen, Prefix = Prefix)
Int = Sys.Int

## Loading the model forcefield into the system
ffstring = parseFF(ffield_file)
Sys.ForceField.SetParamString(ffstring)

## Adding Rg measure
Rg_measure = sim.measure.Rg(Sys, StepFreq = NStep['StepFreq'])
Sys.Measures.extend([Rg_measure])


## Adding a visualizer - only for debugging purposes
if useVis:
	import sim.system.visualize
	Vis = sim.system.visualize.Visualizer3D(Sys, ShowBonds = True)
	Vis.AddAction(Int, TimeFreq = 1.0)

# Initialization
sim.system.init.positions.CubicLattice(Sys, Random = 0.1)
sim.system.init.velocities.Canonical(Sys, Temp = TempSet)


###############     DIFFERENT MD STEPS    ######################
	
## 1. Equilibration
print "Equilibrating system with Langevin Thermostat..."
Int.Method = Int.Methods.VVIntegrate
Int.Method.TimeStep = NStep['Stepsize']

# LD-only case only works with the Nose-Hoover thermostat,
# Integration fails for the Langevin theromostat even with
# very low values of LangevinGamma
if cgtype == 'LD':
	Int.Method.Thermostat = Int.Method.ThermostatNoseHoover
else:	
	Int.Method.Thermostat = Int.Method.ThermostatLangevin
	Int.Method.LangevinGamma = 0.01

Sys.Measures.Reset()
Sys.Measures.VerboseOutput(StepFreq = NStep["StepFreq"])
Sys.TempSet = TempSet
Int.Reset()
Int.Run(NStep["Equil"], 'Equilibrating')


## 2. Production
print "Production..."
Int.Method.Thermostat = Int.Method.ThermostatNoseHoover
Sys.Measures.Reset()
Sys.Measures.VerboseOutput(StepFreq = NStep["StepFreq"])

# Add a data logging action
PElogfile = Prefix + ".log"    
log_action = writeToFile(PElogfile, Sys = Sys)                         
Int.AddAction(log_action)
    
## Add a Lammps format trajectory dump action
Trj = sim.traj.lammps.LammpsWrite(FileName = Prefix + '.lammpstrj.gz', MinImage = True)
Trj.AddAction(Int, StepFreq = NStep['StepFreq'])

Int.Reset()  
Int.Run(NStep["Prod"], 'Production runs') 

## Cleanup
Int.DelAction(log_action)
Trj.DelAction(Int)
Trj.Close()

#!/usr/bin/env python

import numpy as np
import os
import sim
import sys


sys.path.append(os.path.expanduser('~/c25ld'))
import makesys
from chem_data import getModelData

d = getModelData()
poly_dict = d['poly_dict']
cut_dict = d['cut_dict']
TempSet = d['TempSet']

## Iteration numbers
NBlocks = 4
NStep = {"Stepsize": 0.003, "Relax": 100000, "StepFreq": 100, "Equil": 1000000, "Prod": NBlocks*5000000}

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
    
    

def writeToFile(PElog, UWlog, Sys):
    """ 
    Integrator action that logs total potential
    and umbrella run details
    """		
    def PrintHead(Sys):
        Sys.Measures.Reset()
        with open(PElog, "w") as fPE:
        	fPE.write("%7s " % "Step" + Sys.Measures.ReportHead()+"\n")
       	with open(UWlog, "w") as fUW:
	       	fUW.write("%s\t%s\t%s\t%s\n" % ("# Step","Rg", "E_harmonic", "x0_Rg"))

    def PrintMeasures(Sys):
        with open(PElog, "a") as fPE:
        	fPE.write("%7d " % Sys.Int.StepNum + Sys.Measures.ReportVals() + "\n")
        with open(UWlog, "a") as fUW:
		this_Rg = Sys.Measures[-1].Val
		Sys.E_harmonic = Sys.ForceField[-1].Val(this_Rg)
        	fUW.write("%d\t%f\t%f\t%f\n" % (Sys.Int.StepNum, this_Rg, Sys.E_harmonic, Sys.Rgcenter))


    action = sim.integrate.Action(StepFreq = NStep["StepFreq"], Fn = PrintMeasures, InitFn = PrintHead)
    return action


## Input from user
## Input format:  python cg_umbrella_MD.py ffield_file fftype Prefix Rg_0 forceconst
ffield_file = sys.argv[1]
Prefix = sys.argv[2]
fftype = sys.argv[3]
boxlen = sys.argv[4]
Rg_0 = float(sys.argv[5])
FConst = float(sys.argv[6])
cgtype = Prefix.split('_')[-2]

## Setting parameters taken from input
Bias = {'FConst': FConst, 'Rg_0': Rg_0}

## Generating *non* energy minimized system
Sys = makesys.MakeSys(fftype, BoxL = boxlen, Prefix = Prefix, isUmbrella = True)
Int = Sys.Int

## Loading the model forcefield into the system
ffstring = parseFF(ffield_file)
Sys.ForceField.SetParamString(ffstring)

## Adding the biasing forcefield and Rg measurement (only *after* loading the CG forcefield)
print 'Adding Biasing potential with Bias: ', Bias
atomtype = sim.chem.AtomType("P", Mass = poly_dict["Mass_mon"], Color = (1,0,0), Radius = 0.5)
RgFilter = sim.atomselect.PolyFilter([atomtype])
Pumbrella = sim.potential.RgUmbrella(Sys, Filter = RgFilter, FConst = Bias['FConst'], RgSet = Bias['Rg_0'], Label = "U")
Sys.ForceField.extend([Pumbrella])
Sys.__setattr__("Rgcenter", Bias['Rg_0'])
Sys.__setattr__("E_harmonic", 0.)
Rg_measure = sim.measure.Rg(Sys, StepFreq = NStep['StepFreq'])

Sys.Measures.extend([Rg_measure])
Sys.ForceField.extend([Pumbrella])


## Compile and energy minimze system
print 'Post-compilation of system object...'
Sys.Load()
sim.system.init.positions.CubicLattice(Sys, Random = 0.1)
sim.system.init.velocities.Canonical(Sys, Temp = TempSet)
for Method in Int.Methods:
	if hasattr(Method, "TimeStep"): Method.TimeStep = 3.e-3
Int.Method = Int.Methods.VVQuench
Int.Run(NStep['Relax'], "Minimizing")


###############     DIFFERENT MD STEPS    ######################

## 1. Equilibration
print "Equilibrating system..."
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
Int.Run(NStep["Equil"])


## 2. Production (Umbrella sampling starts here)
print "Production..."
Int.Method.Thermostat = Int.Method.ThermostatNoseHoover
Sys.Measures.Reset()
Sys.Measures.VerboseOutput(StepFreq = NStep["StepFreq"])

#Add a data-logging action
PElogfile = Prefix + ".log"
UWlogfile = Prefix + ".colvars.traj"    
log_action = writeToFile(PElogfile, UWlogfile, Sys = Sys)                         
Int.AddAction(log_action)    

## Add a Lammps format trajectory dump action
Trj = sim.traj.lammps.LammpsWrite(FileName = Prefix + 'lammpstrj.gz', MinImage = True)
Trj.AddAction(Int, StepFreq = NStep['StepFreq'])

Int.Reset()  
Int.Run(NStep["Prod"]) 

## Cleanup
Int.DelAction(log_action)
Trj.DelAction(Int)
Trj.Close()
    

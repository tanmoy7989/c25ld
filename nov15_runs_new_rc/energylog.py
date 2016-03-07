#!/usr/bin/env python

import numpy as np
import os, sys, pickle, copy
import sim

#dependencies
import pickleTraj
import parse_potential as pp
sys.path.append('/home/cask0/home/tsanyal/c25ld/srel_rc')
import alg

DEBUG = False
calcAtomEne = True

# user input
fftype = sys.argv[1]
isPolymer = bool(int(sys.argv[2]))
hasLJ = bool(int(sys.argv[3]))
N_mon = int(sys.argv[4])
N_water = int(sys.argv[5])
LammpsTraj = sys.argv[6] 
Prefix = sys.argv[7]
cgffPrefix = sys.argv[8]
LDCut = float(sys.argv[9])

if calcAtomEne: Prefix += '_atom'
    

#check if already computed
LogFileName = '/home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/geom_prop/energy_logs/%s_energy.pickle' % Prefix
#if os.path.isfile(LogFileName): exit()

#system parameters
SysData = alg.makeSysData(isPolymer = isPolymer, hasLJ = hasLJ, N_mon = N_mon, N_water = N_water)
TempSet = 298.0
Delta = 1.2

#CG forcefield parameters
cg_ff_file = cgffPrefix + '_sum.txt'
paramstring = pp.parseParamString(cg_ff_file)

#boxlength extracted from AA traj
Trj = pickleTraj(LammpsTraj)
BoxL = Trj.FrameData['BoxL']

#create the system object
Sys = alg.makeSys(LDCut = LDCut, SysData = SysData, fftype = fftype, BoxL = BoxL, Prefix = Prefix, paramstring = paramstring, Delta = 1.2)
Sys.TempSet = TempSet
Int = Sys.Int
Int.Method = Int.Methods.VVIntegrate
Int.Method.Thermostat = Int.Method.ThermostatLangevin
Int.Method.LangevinGamma = 0.01

#if calcAtomEne and DEBUG: Sys_test = copy.deepcopy(Sys)

# Extract potential objects
Plj = None; Pspline = None ; Plocaldensity = None
if hasLJ:
    Plj = Sys.ForceField[0]
    Pspline = Sys.ForceField[1]    
else:
    Pspline = Sys.ForceField[0]
Plocaldensity = Sys.ForceField[-1]
Pangle = Sys.ForceField[-2]
Pbond = Sys.ForceField[-3]

# Artifically parameterize the LJ potential
Pbond.SetParam(FConst = 0., Dist0 = 0.)
Pangle.SetParam(FConst = 0., Theta0 = 0.)
Plj.SetParam(Epsilon = SysData['poly_dict']['LJEpsilon'], Sigma = SysData['poly_dict']['LJSigma'])
if fftype == 'lj': LJCut = SysData['cut_dict']['LJCut'] 
else: LJCut = SysData['cut_dict']['WCACut']
Plj.Cut = LJCut

# loop over traj frames
start = 0; stop = len(Trj); freq = 1
FrameRange = range(start, stop, freq)
Nframes = len(FrameRange)
Ene_Ind = 0
EnePairSpline = np.zeros(Nframes)
EnePairLJ = np.zeros(Nframes)
EnePairLD = np.zeros(Nframes)
EneSolv_atom = np.zeros([Nframes, N_mon])

pb = sim.utility.ProgressBar(Text = 'Processing frames...',  Steps = Nframes) 
Sys.Flags.Alg.CalcTerms = True   
for frame in FrameRange:
    Pos = Trj[frame]
    Sys.Arrays.Pos = Pos
    
    if not calcAtomEne:
        Sys.ForceField.Eval()
        if DEBUG: 
            print 'Frame = %d, PairSpline = %11.4e, PairLocalDensity = %11.4e, LJ = %11.4e\n' % (Ene_Ind, EnePairSpline[Ene_Ind], EnePairLD[Ene_Ind], EnePairLJ[Ene_Ind])
            raw_input()
        Ene = Sys.ForceField.Terms
        EnePairLJ[Ene_Ind] = Ene[0]
        EnePairSpline[Ene_Ind] = Ene[1]
        EnePairLD[Ene_Ind] = Ene[-1]
   
    else:
        if DEBUG: 
            part1 = np.arange(2); part2 = np.arange(2,25)
            Sys.ForceField.Eval(part1); x1 = copy.deepcopy(Sys.ForceField.Terms)
            Sys.ForceField.Eval(part2); x2 = copy.deepcopy(Sys.ForceField.Terms)
            Sys.ForceField.Eval(); x_all = copy.deepcopy(Sys.ForceField.Terms)
            print  "Frame = ", Ene_Ind
            print  "Part1 = ", x1
            print  "Part2 = ", x2
            print  "Sum12 = ", x1+x2
            print  "Whole = ", x_all
            print "\n-----------------------------------\n\n"
            raw_input()
        
        else:
            for n in range(N_mon):
                Sys.ForceField.Eval(np.arange(n+1))
                Ene = Sys.ForceField.Terms
                EneSolv_atom[Ene_Ind, n] = (Ene[1] - Ene[0]) + Ene[-1]
    
    Ene_Ind += 1
    pb.Update(Ene_Ind)

# dump data
if not calcAtomEne: data = (EnePairSpline, EnePairLD, EnePairLJ)
else: data = EneSolv_atom.flatten()
pickle.dump(data, open(LogFileName, 'w'))
        

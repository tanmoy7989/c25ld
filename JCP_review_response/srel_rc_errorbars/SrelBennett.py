#!/usr/bin/env python

import os, sys, pickle
import numpy as np
import sim
import pickleTraj
import parse_potential

#import subroutine module
sys.path.append('/home/cask0/home/tsanyal/c25ld/srel_rc')
import alg

#flags
DEBUG = False
RefreshLogs = True
doBlockAvg = True
NBlocks = 10
MeasureFreq = 1

def EvalSysEne(Sys, sampleTrj, serial_index, TempSet = 298.0, Iter = (0,-1,100)):
    '''
    sampleTrj refer to the state (AA or CG) from which the position
    co-ordintates are being sampled. Sys is a system with the forcefield
    of the state at which energies are to be evaluated
    returns beta * Penergy of system
    '''
    
    start = Iter[0]; stop = Iter[1]; freq = Iter[2]
    if stop == -1: stop = len(sampleTrj)  
    ind = serial_index[Sys.sample_trajtype]
    
    FrameRange = range(start, stop, freq)
    Nframes = len(FrameRange)
    Ene = np.zeros([Nframes], np.float64)
    beta = 1./(TempSet * Sys.Units.kB)
    
    #for each frame, minimage and calculate energy
    Ene_Ind = 0
    pb = sim.utility.ProgressBar(Text = 'Processing frames...',  Steps = Nframes)
    for frame in FrameRange:
        Pos = sampleTrj[frame]
        Pos = Pos[ind[0]:ind[1], :]
        
        Sys.Arrays.Pos = Pos
        Sys.ForceField.Eval()
        Ene[Ene_Ind] = Sys.PEnergy
  
        if DEBUG:
            print frame, Ene_Ind, Ene[Ene_Ind]
            raw_input()
            
        pb.Update(Ene_Ind)
        Ene_Ind += 1

    if DEBUG: 
        print '(sample, eval)', Sys.sample_trajtype, Sys.eval_trajtype
        print len(Ene), len(sampleTrj)
        raw_input()
    
    return beta * Ene
        
   
           
####### MAIN #######

'''USAGE: 
sample command line call 
python dSreldrc.py lj 0 0 25 1700 /home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane/nvt/methane_lj_unbiased_mapped.lammpstrj.gz methane_%s_%d_SPLD_sum.txt methane_%s_%d_MD.lammpstrj.gz methane_lj 
'''

#user input
fftype = sys.argv[1]
isPolymer = bool(int(sys.argv[2]))
hasLJ = bool(int(sys.argv[3]))
N_mon = int(sys.argv[4])
N_water = int(sys.argv[5])
LammpsTraj = sys.argv[6] 
ff_fmtstring = sys.argv[7]
md_fmtstring = sys.argv[8]
datalog_prefix = sys.argv[9]

#system parameters
SysData = alg.makeSysData(isPolymer = isPolymer, hasLJ = hasLJ, N_mon = N_mon, N_water = N_water)
TempSet = 298.0
Delta = 1.2

#input files
base_fmt = {'cg_ff': ff_fmtstring, 'cg_trj': md_fmtstring}

#output arrays
LDCuts = np.arange(4.5, 10.0, 0.1)
delta_Srel = {'Bennett': np.zeros([len(LDCuts), NBlocks], np.float64), 'Fep' : np.zeros([len(LDCuts), NBlocks], np.float64)}
datalog = {'Bennett': '%s_Srel_Bennett.dat' % datalog_prefix, 'Fep': '%s_Srel_Fep.dat' % datalog_prefix }

#various indexing schemes
poly_dict = SysData['poly_dict']
serial_index = {'AA': (0, poly_dict['N_mon']*poly_dict['N_poly']), 
                'CG': (0, poly_dict['N_mon']*poly_dict['N_poly'])} #always using a mapped traj

#extract all necessary AA data
Trj_AA = pickleTraj(LammpsTraj)
BoxL = Trj_AA.FrameData['BoxL']

#begin looping over different LD cutoffs
for i, LDCut in enumerate(LDCuts[0:-1]):
    
    print '\n\nLDCUTs = (%g, %g)\n' % (LDCuts[i], LDCuts[i+1])
  
    #extract all necessary CG data
    sumfile1 = base_fmt['cg_ff'] % (fftype, i); sumfile2 = base_fmt['cg_ff'] % (fftype, (i+1))
    trajfile1 = base_fmt['cg_trj'] % (fftype, i) ; trajfile2 = base_fmt['cg_trj'] % (fftype, (i+1))
    ParamString1 = parse_potential.parseParamString(sumfile1); ParamString2 = parse_potential.parseParamString(sumfile2)
    Trj_CG_1 = pickleTraj(trajfile1); Trj_CG_2 = pickleTraj(trajfile2)
    
    #create systems with forcefields of LDCuts i and i+1
    Sys1 = alg.makeSys(LDCut = LDCuts[i], Delta = Delta, BoxL = BoxL, fftype = fftype, Prefix = 'bennettsys1',
                       paramstring = ParamString1, sample_trajtype = 'CG', eval_trajtype = 'CG', SysData = SysData)
                       
    Sys2 = alg.makeSys(LDCut = LDCuts[i+1], Delta = Delta, BoxL = BoxL, fftype = fftype, Prefix = 'bennettsys2',
                       paramstring = ParamString2, sample_trajtype = 'CG', eval_trajtype = 'CG', SysData = SysData)
                      
    # block average
    if not doBlockAvg: NBlocks = 1
    for j, b in enumerate(range(NBlocks)):
    	print '\n======== BLOCK %d ==============\n' % b
    	
    	block_start_AA = b * int(len(Trj_AA)/NBlocks)
    	block_end_AA = (b+1) * int(len(Trj_AA)/NBlocks)
    	
    	block_start_CG = b * int(len(Trj_CG_1)/NBlocks)
    	block_end_CG = (b+1) * int(len(Trj_CG_1)/NBlocks)

    	TrjIter = {'AA': [block_start_AA, block_end_AA, MeasureFreq], 'CG': [block_start_CG, block_end_CG, MeasureFreq]}

    	#calculate Energy difference in AA ensemble
    	print '--> Calculating energies in AA ensemble'
    	Sys1.sample_trajtype = 'AA'; Sys2.sample_trajtype = 'AA'
    	BE1_AA = EvalSysEne(Sys = Sys1, sampleTrj = Trj_AA, Iter = TrjIter['AA'], serial_index = serial_index)
    	BE2_AA = EvalSysEne(Sys = Sys2, sampleTrj = Trj_AA, Iter = TrjIter['AA'], serial_index = serial_index)
    
    	#reprocess trajectories to calculate free energy differences
    	print '--> Reprocessing trajectories'
    	Sys1.sample_trajtype = 'CG'; Sys2.sample_trajtype = 'CG'
    
    	BE11 = EvalSysEne(Sys = Sys1, sampleTrj = Trj_CG_1, Iter = TrjIter['CG'], serial_index = serial_index)
    	BE21 = EvalSysEne(Sys = Sys1, sampleTrj = Trj_CG_2, Iter = TrjIter['CG'], serial_index = serial_index)
    	BE22 = EvalSysEne(Sys = Sys2, sampleTrj = Trj_CG_2, Iter = TrjIter['CG'], serial_index = serial_index)
    	BE12 = EvalSysEne(Sys = Sys2, sampleTrj = Trj_CG_1, Iter = TrjIter['CG'], serial_index = serial_index)
    
    	Nframes1 = len(BE11); Nframes2 = len(BE22)
    
    	#running Bennett's algorithm
    	print "--> Running Bennett's method..."
    	FE1, Err1 = alg.BennettFE(BE11 = BE11, BE22 = BE22, BE12 = BE12, BE21 = BE21, Verbose = False)
    	delta_Srel1 = (np.mean(BE2_AA) - np.mean(BE1_AA)) - FE1
   
    	print "--> Running 1D Free Energy Perturbation..."
    	FE2 = alg.FEP(BE11 = BE11, BE22 = BE22, BE12 = BE12, BE21 = BE21, Verbose = False)
    	delta_Srel2 = (np.mean(BE2_AA) - np.mean(BE1_AA)) - FE2
    	Err2 = 0.0
    	
    	for key in datalog.keys():
    	   if key == 'Bennett':
    	       delta_Srel[key][i,j] = delta_Srel1
    	   elif key == 'Fep':
    	       delta_Srel[key][i,j] = delta_Srel2


# write to file
for key in datalog.keys():
	np.savetxt(datalog[key], delta_Srel[key])

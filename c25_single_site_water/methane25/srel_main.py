#!/usr/bin/env python

import os, sys, pickle
import parse_potential as pp

sys.path.append('/home/cask0/home/tsanyal/c25ld/c25_single_site_water')
import cgmodel as cg

fftype = 'wca'
isPolymer = 0
LammpsTraj = '/home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane25/nvt/methane_wca_unbiased_single_site_water.lammpstrj.gz'
N_mon = 25
N_water = 1700
LDCut = 7.8 #based on 2nd shell cutoff from monomer-water rdf
Prefix = 'methane25_wca_ss'

SysData = cg.makeSysData(isPolymer = isPolymer, N_mon = N_mon, N_water = N_water)
Sys = cg.makeSys(LDCut = LDCut, SysData = SysData, Prefix = Prefix)

### SREL CG ###
opt_cases = ['SP', 'SPLD']
cg.runSrel(Sys = Sys, SysData = SysData, LammpsTraj = LammpsTraj, fixWW = True)

### CG MD ###
#for x in opt_cases:
#        if not os.path.isfile(Prefix + '_%s.lammpstrj.gz' % x):
#            s = pp.parseParamString(Prefix + '_%s_sum.txt' % x)
#            cg.runMD(Sys = Sys, SysData = SysData, ParamString = s, Prefix = Prefix + '_%s' % x )
